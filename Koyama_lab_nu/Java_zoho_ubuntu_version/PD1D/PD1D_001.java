//*** [プログラム (PD1D_001.java)] ************************************************
//*** [インポート文] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class PD1D_001 extends Frame{

//*** [グローバル変数] *****************************************************************************************
	static int ND=512;		//組織１辺の分割数
	static int nd=ND;			//濃度の分割数
	static int ndm=ND-1;	//濃度の分割数-1
	static int width;			// Window全体の幅
	static int height;			// Window全体の高さ
	static int xwidth;			// 描画領域の幅
	static int yheight;			// 描画領域の高さ
	static int insetx;			// Windowの枠の幅（左右および下）
	static int insety;			// Windowの枠の幅（上）
	static double PI=3.141592;				//π
	static double RR=8.3145;				//ガス定数
	static double [] ch=new double[ND];	//組織内の濃度デ－タ配列
	static Graphics g;						//自由エネルギー曲線画面のグラフィックスオブジェクト
	static double time1;						//計算時間（カウント数）
	static double temp; 						//温度(K)
	static double c0;							//合金組成（モル分率）

//*** [コンストラクタ] ****************************
	public PD1D_001(){
		xwidth=800; yheight=200; 			//描画画面の横と縦の長さ（ピクセル単位）
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

//*** [メインプログラム] *******************************************************************************************
	public static void main(String[] args) throws Exception{//例外処理は行わない

		PD1D_001 prog=new PD1D_001();	//PD1D_001のインスタンスprogを生成

		int i;										//整数
		int ip, im; 								//整数（i+1, i-1）
		double delt;								//時間刻み（無次元）
		double al;									//計算領域の１辺の長さ
		double time1max;						//計算時間（最大カウント数）
		double [] ck = new double[ND];		//拡散ポテンシャル
		double [] ch2 = new double[ND];		//組織内の濃度デ－タ予備配列
		double mu_chem, mu_surf;			//各ポテンシャル
		double c1, c2;							//濃度(A:1, B:2)
		double L0;								//原子間相互作用パラメータ
		double kappa_c;							//濃度勾配エネルギー係数
		double Mc;								//易動度
		double c_flu;								//濃度場の揺らぎの大きさ
		double cddtt;								//濃度の増分

		double b1;								//差分ブロックサイズ
		double c2ip, c2im; 						//差分においてc2を中心に、その左右の濃度
		double sumc, dc; 						//濃度場の総和、平均組成からのずれ

//---- 各種パラメータ設定 ----------------------------------------------------------------------
		temp=1000.0;			//[K]
		c0=0.1;					//合金の平均組成（ここではA-40at%B合金を設定している）
		delt=0.01;
		time1=0.0;				//計算の繰返し回数
		time1max=1.0e+08;		//繰返し回数の最大値

		al=500.0;					//２次元計算領域の１辺の長さ(nm)
		al=al*1.0e-9;				//(m)に変換
		b1=al/(double)nd;		//差分１ブロックのサイズ

		L0=2.5e+04;				//原子間相互作用パラメータ（J/mol）
		L0=L0/RR/temp;			//無次元化

		kappa_c=5.0e-15;		//濃度勾配エネルギ－係数、単位は[J m^2/mol]
		kappa_c=kappa_c/b1/b1/RR/temp;		//(b1^2*rr*temp)で無次元化

		Mc=c0*(1.0-c0);			//拡散の易動度
		c_flu=0.1;					//濃度ゆらぎ振幅の最大値

//---- 時間0における初期濃度プロファイル設定 ------------------------------------------------------
		prog.ini_comp_field();	//乱数によって生成する場合
		//prog.datin();				//ファイルから読み込む場合

//---- 濃度プロファイルの時間発展の計算 ---------------------------------------------------------------
		while(time1<=time1max){

//---- 濃度プロファイルの表示 -------------------------------------------------------------------------------
			if((((int)(time1) % 500)==0)){ prog.repaint(); }	
											//カウント数が500の倍数おきに濃度プロファイル描画

//---- 濃度プロファイルの保存 ------------------------------------------------------------------------------
			if((((int)(time1) % 1000)==0)){ prog.datsave(); }	
											//カウント数が1000の倍数おきに濃度プロファイル保存
			//if(time1==3000.0){ prog.datsave(); }	//カウント数が3000の時に濃度プロファイル保存

//---- 拡散ポテンシャルの計算 -------------------------------------------------------------------------------
			for(i=0;i<=ndm;i++){
				ip=i+1; im=i-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}						//周期的境界条件
				c2=ch[i];		c1=1.0-c2;									//位置iにおけるc1とc2
				c2ip=ch[ip]; c2im=ch[im];									//差分においてc2の左右の濃度

	 			mu_chem=L0*(c1-c2)+Math.log(c2)-Math.log(c1);		//化学ポテンシャル差
				mu_surf=-2.0*kappa_c*(c2ip+c2im-2.0*c2);				//濃度勾配のポテンシャル
				ck[i]=mu_chem+mu_surf; 									//拡散ポテンシャル
			}

//---- 濃度場の時間変化(非線形拡散方程式の差分解陽解法) -------------------------------------------
			for(i=0;i<=ndm;i++){
				ip=i+1; im=i-1;  if(i==ndm){ip=0;}	if(i==0){im=ndm;} 	//周期的境界条件
				cddtt=Mc*( ck[ip]+ ck[im]-2.0* ck[i]);						//非線形拡散方程式
				//ch2[i]=ch[i]+cddtt*delt; 									//濃度場の時間発展
				ch2[i]=ch[i]+( cddtt+c_flu*(2.0*Math.random()-1.0) )*delt;	//濃度場の時間発展（濃度揺らぎを考慮）
			}

//*** [濃度場の収支の補正] *******************************************************
//*** 数値計算であるので濃度場の収支の補正を行う（実際には毎ステップ行う必要はない）****
  			sumc=0.; for(i=0;i<=ndm;i++){ sumc+=ch2[i]; }				//濃度プロファイルの積分
			dc=sumc/(double)nd-c0;										//濃度プロファイルの変動量

			for(i=0;i<=ndm;i++){
				ch[i]=ch2[i]-dc;												//濃度場の補正
				if(ch[i]>=1.){ch[i]=1.0-1.0e-6;}								//濃度が1を超えた場合の補正
				if(ch[i]<=0.){ch[i]=1.0e-6;}									//濃度が0を切った場合の補正
			}

//******[時間増加]**************************************************************
			time1=time1+1.0;												//計算時間の増加
		}//while

//----------------------------------------------------------------

		System.out.printf("\n 終了しました。グラフの右上×をクリックして終了してください。\n");

}//main

// 以下はサブルーチン
// **** 初期濃度場設定 *****************************************************
	public void ini_comp_field(){
		int i;
		double rnd0, fac1;

		fac1=0.01; 				//初期濃度揺らぎの最大変化量を1%に設定
		for(i=0;i<=ndm;i++){
			rnd0=2.0*Math.random()-1.0;  ch[i]=c0+rnd0*fac1;	//乱数にて初期濃度場を設定
		}
	}

// **** 濃度場描画 *****************************************************
	public void paint(Graphics g){
		int i, ii, i1, ii1, i2, ii2;
		double x, xmax, xmin, dx, y, ymax, ymin, dy, d0;
		double c, cmax, cmin;
		double gx1, gy1, gx2, gy2;
		int ixmin, iymin, ixmax, iymax, igx1, igy1, igx2, igy2, id0;
		int idx, idy;
		double col;

		xmin=0.; xmax=1.; dx=0.1; 								//横軸の最小値、最大値、分割間隔（実際の値）
		ymin=0.; ymax=1.; dy=0.1;									//縦軸の最小値、最大値、分割間隔（実際の値）
		cmin=0.0; cmax=1.0; 										//濃度の最小値、最大値（実際の値）
		ixmin=0; iymin=0;  ixmax=xwidth; iymax=yheight;		//ピクセル単位の場合
		idx=(int)(0.1*ixmax);
		idy=(int)(0.1*iymax);

		g.setColor(Color.white); 									//色を白に設定
		g.fillRect(insetx, insety, xwidth, yheight); 				//画面を上で指定した色で塗る
		g.setColor(Color.black); 									//色を黒に設定
		g.drawRect(insetx, insety, xwidth, yheight); 				//グラフの外周を上で指定した色で描く

		//グラフ内の縦横の線を等間隔に描く
		for(i=0;i<=ixmax;i+=idx){ g.drawLine(insetx+i, insety+iymin, insetx+i, insety+iymax); }
		for(i=0;i<=iymax;i+=idy){ g.drawLine(insetx+ixmin, insety+i, insetx+ixmax, insety+i); }

		d0=1.0/(double)nd/2.0;										//差分ブロックの長さの半分
		id0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*d0 );		//d0のピクセル化

		System.out.printf("%f \n", time1);							//計算の繰返し回数を標準入出力に表示

		g.setColor(Color.red); 										//色を赤に設定
		for(i=0;i<=ndm;i++){
			//濃度プロファイルの値（実際の値）
			i1=i; i2=i+1;
			gx1=1./(double)nd*(double)i1+d0;		gx2=1./(double)nd*(double)i2+d0;
			ii1=i1; ii2=i2;  if(i==ndm){ii2=0;}
			gy1=(ch[ii1]-cmin)/(cmax-cmin); 		gy2=(ch[ii2]-cmin)/(cmax-cmin);

			//濃度プロファイルの値をスクリーン座標に変換
			igx1=(int)( ((double)ixmax-(double)ixmin)*(gx1-xmin)/(xmax-xmin)+(double)ixmin );
			igy1=(int)( (double)iymin+(double)iymax
									-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy1-ymin)+(double)iymin) );
			igx2=(int)( ((double)ixmax-(double)ixmin)*(gx2-xmin)/(xmax-xmin)+(double)ixmin );
			igy2=(int)( (double)iymin+(double)iymax
									-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy2-ymin)+(double)iymin) );

			g.drawLine(insetx+igx1, insety+igy1, insetx+igx2, insety+igy2);	//濃度プロファイルの描画
		}
	}

//*** [デ－タの保存] ************************************
	private void datsave() throws Exception{
		int i;

		//保存ファイル名をtest.datとする。
		PrintWriter outfile= new PrintWriter(
			new BufferedWriter(new FileWriter("test.dat", true)) );	//ファイルのオープン(追記)

		outfile.println(time1);								//カウントの書き込み
		for(i=0;i<=ndm;i++){ outfile.println(ch[i]); }		//濃度場の書き込み
		outfile.close();										//ファイルのクローズ
	}

//*** [デ－タの読込み] ************************************
	private void datin() throws Exception{
		int i;
		String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//ファイルのオープン

		s_data=infile.readLine();							//文字列として読み込み
		time1=new Double(s_data).doubleValue();		//文字を数値へ変換
		for(i=0;i<=ndm;i++){
			s_data=infile.readLine();						//文字列をして読み込み
			ch[i]=new Double(s_data).doubleValue();	//文字を数値へ変換
		}
		infile.close();											//ファイルのクローズ
	}

//****************************************************************
}//PD1D_001
//*** プログラム終了 ************************************************************
