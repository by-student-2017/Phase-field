//*** [プログラム (PD2D_001.java)] ************************************************
//*** [インポート文] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class PD2D_001 extends Frame{

//*** [グローバル変数] *****************************************************************************************
	static int ND=64;			//組織１辺の分割数
	static int nd=ND;			//濃度の分割数
	static int ndm=ND-1;	//濃度の分割数-1
	static int width;			// Window全体の幅
	static int height;			// Window全体の高さ
	static int xwidth;			// 描画領域の幅
	static int yheight;			// 描画領域の高さ
	static int insetx;			// Windowの枠の幅（左右および下）
	static int insety;			// Windowの枠の幅（上）
	static double PI=3.141592;						//π
	static double RR=8.3145;						//ガス定数
	static double [][] ch=new double[ND][ND];	//組織内の濃度デ－タ配列
	static Graphics g;								//自由エネルギー曲線画面のグラフィックスオブジェクト
	static double time1;								//計算時間（カウント数）
	static double temp; 								//温度(K)
	static double c0;									//合金組成（モル分率）


//*** [コンストラクタ] ****************************
	public PD2D_001(){
		xwidth=400; yheight=400; 			//描画画面の横と縦の長さ（ピクセル単位）
		insetx=4; insety=30;					//描画画面の枠の長さ
		width=xwidth+insetx*2;  			//描画Window全体の横の長さ
		height=yheight+insetx+insety; 	//描画Window全体の縦の長さ
		setSize(width, height);				//描画Windowのセット
		setBackground(Color.white); 		//描画Windowの描画部の色を白に設定
		setVisible(true); 						//描画Windowを見えるようにする
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
												//Windowを閉じる時の操作（Windowの右上×の設定）
		});
	}

//*** [main プログラム] *****************************
	public static void main(String[] args) throws Exception{//例外処理は行わない

		PD2D_001 prog=new PD2D_001();	//PD2D_001のインスタンスprogを生成

		int i, j; 										//整数
		int ip, im, jp, jm; 								//整数（i+1, i-1, j+1, j-1）
		double delt;									//時間刻み（無次元）
		double al;										//計算領域の１辺の長さ
		double time1max;							//計算時間（最大カウント数）
		double [][] ck = new double[ND][ND];		//拡散ポテンシャル
		double [][] ch2 = new double[ND][ND]; 	//組織内の濃度デ－タ予備配列
		double mu_chem, mu_surf;				//各ポテンシャル
		double c1, c2;								//濃度(Fe:1, Cu:2)
		double L0;									//原子間相互作用パラメータ
		double kappa_c;								//濃度勾配エネルギー係数
		double Mc;									//易動度関数とその微分
		double c_flu;									//濃度場の揺らぎの大きさ
		double cddtt;									//濃度の増分

		double b1;									//差分ブロックサイズ
		double c2ip, c2im, c2jp, c2jm; 			//差分ブロックにおいてc2を中心に、その上下左右の濃度
		double sumc, dc; 							//濃度場の総和、平均組成からのずれ


//---- 各種パラメータ設定 ----------------------------------------------------
		temp=1000.0;  			//[K]
		c0=0.4;					//合金の平均組成（ここではA-30at%B合金を設定している）
		delt=0.04;					//時間きざみ
		time1=0.0;				//計算の繰返し回数
		time1max=1.0e+08;		//繰返し回数の最大値

		al=60.0; 					//２次元計算領域の１辺の長さ(nm)
		al=al*1.0e-9;				//(m)に変換
		b1=al/(double)nd;		//差分１ブロックのサイズ

		L0=2.5e+04;				//原子間相互作用パラメータ（J/mol）
		L0=L0/RR/temp;			//無次元化

		kappa_c=5.0e-15;			//濃度勾配エネルギ－係数、単位は[J m^2/mol]
		kappa_c=kappa_c/b1/b1/RR/temp;		//濃度勾配エネルギ－係数、(b1^2*rr*temp)で無次元化

		Mc=c0*(1.0-c0);			//拡散の易動度
		c_flu=0.1;

//---- 時間0における初期濃度場設定 ----------------------------------------------------
		prog.ini_comp_field();
		//prog.datin();				//ファイルから読み込む場合

//---- 濃度場の時間発展の計算 ----------------------------------------------------
		while(time1<=time1max){

//---- 濃度場の表示 ----------------------------------------------------
			//カウント数が200の倍数おきに濃度場描画
			if((((int)(time1) % 200)==0)){ prog.update_draw(g); }	//描画時のチラツキを抑えるため
			//if((((int)(time1) % 200)==0)){ prog.repaint(); }

//---- 濃度場の保存 ----------------------------------------------------
			if((((int)(time1) % 500)==0)){ prog.datsave(); }			//カウント数が500の倍数おきに濃度場保存
			//if(time1==3000.0){ prog.datsave(); }						//カウント数が3000の時に濃度場保存

//---- 拡散ポテンシャルの計算 -----------------------------------------------------
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1; jp=j+1; jm=j-1;
					if(i==ndm){ip=0;}  if(i==0){im=ndm;} 	//周期的境界条件
					if(j==ndm){jp=0;}  if(j==0){jm=ndm;} 	//周期的境界条件

					c2=ch[i][j]; 		c1=1.0-c2; 				//位置(i,j)におけるc1とc2
					c2ip=ch[ip][j]; c2im=ch[im][j]; c2jp=ch[i][jp]; c2jm=ch[i][jm]; 	//差分においてc2の前後左右の濃度

		 			mu_chem=L0*(c1-c2)+Math.log(c2)-Math.log(c1); 				//化学ポテンシャル差
					mu_surf=-2.0*kappa_c*(c2ip+c2im+c2jp+c2jm-4.0*c2);		//濃度勾配のポテンシャル
					ck[i][j]=mu_chem+mu_surf; 										//拡散ポテンシャル
				}
			}

//---- 濃度場の時間変化(非線形拡散方程式の差分解陽解法) ------------------------------------
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
					if(i==ndm) {ip=0;}	if(i==0) {im=ndm;} 	//周期的境界条件
					if(j==ndm) {jp=0;}	if(j==0) {jm=ndm;} 	//周期的境界条件
					cddtt=Mc*(ck[ip][j]+ck[im][j]+ck[i][jp]+ck[i][jm]-4.0* ck[i][j]);		//非線形拡散方程式
					//ch2[i][j]=ch[i][j]+cddtt*delt; 										//濃度場の時間発展
					ch2[i][j]=ch[i][j]+( cddtt+c_flu*(2.0*Math.random()-1.0) ) *delt;	//濃度場の時間発展（濃度揺らぎ導入）
				}
			}

//*** [濃度場の収支の補正] *******************************************************
//*** 数値計算であるので濃度場の収支の補正を行う（実際には毎ステップ行う必要はない）。] ****
  		sumc=0.;
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					sumc=sumc+ch2[i][j]; 					//濃度場の積分
				}
			}
			dc=sumc/(double)nd/(double)nd-c0;			//濃度場の変動量

			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ch[i][j]=ch2[i][j]-dc; 						//濃度場の補正
					if(ch[i][j]>=1.){ch[i][j]=1.0-1.0e-6;} 	//濃度が1を超えた場合の補正
					if(ch[i][j]<=0.){ch[i][j]=1.0e-6;} 			//濃度が0を切った場合の補正
				}
			}

//******[時間増加]*************************************************
			time1=time1+1.0; 	//計算時間の増加
		}//while

//----------------------------------------------------------------

		System.out.printf("\n 終了しました。グラフの右上×をクリックして終了してください。\n");

}//main

// 以下はサブルーチンである。
// **** 初期濃度場設定 *****************************************************
	public void ini_comp_field(){
		int i, j;
		double rnd0, fac1;

		fac1=0.01; 	//初期濃度揺らぎの最大変化量を1%に設定
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				rnd0=2.0*Math.random()-1.0;  ch[i][j]=c0+rnd0*fac1; 	//乱数にて初期濃度場を設定
			}
		}

	}

// *********************************************************
//チラツキ防止のためdraw関数を別に定義
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** 濃度場描画 *****************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);		//Windowをクリア

		int i, j, ii, jj;
		double c, x, xmax, xmin, y, ymax, ymin, rad0;
		int ixmin=0, iymin=0, igx, igy, irad0;
		int ixmax=xwidth, iymax=yheight;
		int icol;

		xmin=0.; xmax=1.; 						//横軸の最小値、最大値
		ymin=0.; ymax=1.; 						//縦軸の最小値、最大値

		rad0=1.0/(double)nd/2.0; 				//差分ブロックの長さの半分
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 ); 	//rad0のピクセル化

		System.out.printf("%f \n", time1); 	//計算の繰返し回数を標準入出力に表示

		for(i=0;i<=nd;i++){
			for(j=0;j<=nd;j++){
				//濃度場の位置座標（実際の値）
				x=1.0/(double)nd*(double)i+rad0;
				y=1.0/(double)nd*(double)j+rad0;
				//濃度場の位置座標（スクリーン座標に変換）
				igx=(int)( ((double)ixmax-(double)ixmin)*(x-xmin)/(xmax-xmin)+(double)ixmin );
				igy=(int)( ((double)iymax-(double)iymin)*(y-ymin)/(ymax-ymin)+(double)iymin );

				//個々の差分ブロックの濃度値
				ii=i; jj=j;
				if(i==nd){ii=0;} if(j==nd){jj=0;}											//周期的境界条件
				icol=(int)(255.0*(1.0-ch[ii][jj]));											//色の諧調をグレースケールにする
				//icol=(int)(255.0*ch[ii][jj]);												//明暗を反転する場合
				if(icol>=255){icol=255;} if(icol<=0){icol=0;}							//明暗の範囲の補正
				g.setColor(new Color(icol,icol,icol)); 									//濃度を明暗で設定
				g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);		//個々の差分ブロック描画
			}
		}
	}

//*** [デ－タの保存] ************************************
	private void datsave() throws Exception{
		int	i, j;

		//保存ファイル名をtest.datとする。
		PrintWriter outfile= new PrintWriter(
			new BufferedWriter(new FileWriter("test.dat", true)) );		//ファイルのオープン追記

		outfile.println(time1);				//カウントの書き込み
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				outfile.println(ch[i][j]);		//濃度場の書き込み
			}
		}
		outfile.close();						//ファイルのクローズ
	}

//*** [デ－タの読込み] ************************************
	private void datin() throws Exception{
		int	i, j;
		String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//ファイルのオープン

		s_data=infile.readLine();  									//文字列をして読み込み
		time1=new Double(s_data).doubleValue();				//文字を数値へ変換
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				s_data=infile.readLine();  							//文字列をして読み込み
				ch[i][j]=new Double(s_data).doubleValue();		//文字を数値へ変換
			}
		}
		infile.close();													//ファイルのクローズ
	}

//****************************************************************
}//PD2D_001
//*** プログラム終了 ************************************************************
