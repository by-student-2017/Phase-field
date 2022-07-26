//*** [プログラム (MT2D_001.java)] ************************************************
//*** [インポート文] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class MT2D_001 extends Frame{
//*** [グローバル変数] *****************************************************************************************
	static int ND=128;				//組織１辺の分割数
	static int IG=7;					//2^IG=ND
	static int nd=ND;					//計算領域の１辺の分割数
	static int ndm=ND-1;			//計算領域の１辺の分割数-1
	static int nd2=ND/2;				//計算領域の１辺の分割数/2
	static int ig=IG;					//2^ig=ND
	static int width;					// Window全体の幅
	static int height;					// Window全体の高さ
	static int xwidth;					// 描画領域の幅
	static int yheight;					// 描画領域の高さ
	static int insetx;					// Windowの枠の幅（左右および下）
	static int insety;					// Windowの枠の幅（上）
	static double PI=3.141592;		//π
	static double RR=8.3145;		//ガス定数
	static double [][] s1h=new double[ND][ND];		//組織内の秩序変数デ－タs1配列
	static double [][] s2h=new double[ND][ND];		//組織内の秩序変数デ－タs2配列
	static Graphics g;									//自由エネルギー曲線画面のグラフィックスオブジェクト
	static double time1;									//計算時間（カウント数）
	static double temp; 									//温度[K]

	//高速フーリエ変換のプログラムについては、文献(9)をご参照下さい。
	static double qs;										//フ－リエ変換(qs:-1)と逆フ－リエ変換(qs:1)の区別
	static double [][] xr = new double[ND][ND];		//フ－リエ変換の実数パ－トに使用する配列
	static double [][] xi = new double[ND][ND];		//フ－リエ変換の虚数パ－トに使用する配列
	static double [] xrf = new double[ND];				//フ－リエ変換の実数パ－トに使用する配列
	static double [] xif = new double[ND];				//フ－リエ変換の虚数パ－トに使用する配列
	static double [] s = new double[ND];				//sinのテーブル
	static double [] c = new double[ND];				//cosのテーブル
	static int [] ik = new int[ND];						//ビット反転操作の配列

//*** [コンストラクタ] ****************************
	public MT2D_001(){
		xwidth=400; yheight=400; 			//描画画面の横と縦の長さ（ピクセル単位）
		insetx=4; insety=30;					//描画画面のふちの長さ
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

		MT2D_001 prog=new MT2D_001();		//MT2D_001のインスタンスprogを生成

		double [][] ec11 = new double[ND][ND];			//拘束歪変動量
		double [][] ec22 = new double[ND][ND];			//拘束歪変動量
		double [][] ep11h0 = new double[ND][ND];		//変態歪
		double [][] ep22h0 = new double[ND][ND];		//変態歪
		double [][] ep11qrh0 = new double[ND][ND];		//拘束歪変動量のフーリエ変換（実数部）
		double [][] ep11qih0 = new double[ND][ND];		//拘束歪変動量のフーリエ変換（虚数部）
		double [][] ep22qrh0 = new double[ND][ND];		//拘束歪変動量のフーリエ変換（実数部）
		double [][] ep22qih0 = new double[ND][ND];		//拘束歪変動量のフーリエ変換（虚数部）
		double [][] eta_s1 = new double[4][4];				//バリアント１の変態歪
		double [][] eta_s2 = new double[4][4];				//バリアント２の変態歪
		double [][] s1k_su = new double[ND][ND];		//勾配ポテンシャル
		double [][] s2k_su = new double[ND][ND];		//勾配ポテンシャル

		double s1, s2;										//マルテンサイトのphase field
		double s1k_chem, s1k_str;							//化学ポテンシャル、弾性ポテンシャル
		double s2k_chem, s2k_str;							//化学ポテンシャル、弾性ポテンシャル
		double c11, c12, c44, lam0, mu0, nu0;			//弾性定数
		double el_fac;										//弾性定数規格化変数
		double ep11T, ep22T;								//個々の歪成分の和
		double ep11_0, ep22_0;							//組織内の変態歪の平均値
		double ep11_a, ep22_a, ep12_a, ep21_a;		//外力に起因する歪
		double sig11_a, sig22_a;							//外力
		double Z11ep, Z12ep, Z21ep, Z22ep;			//フーリエ変換時に使用する係数
		double sum11, sum22;
		double s1ddtt, s2ddtt;								//phase fieldの時間変化量
		double delt;											//時間きざみ

		int   i, j, k, l, ii=0, jj=0, kk, iii, jjj;								//整数
		int   p, q, m, n;										//整数
		int   ip, im, jp, jm;									//整数
		double al, temp;										//計算領域、温度
		double time1max;									//最大時間（計算を止める際に使用）
		double b1, vm0, atom_n;							//規格化長さ、モル体積、単位胞内の原子数
		double smob;										//マルテンサイト変態ダイナミクスの緩和係数
		double nx, ny, nxx, nyy, alnn;						//逆空間の基本ベクトル、その自乗、ノルム

		double AA0, AA1, AA2, AA3;						//化学的駆動力定数
		double a1_c, b1_c, c1_c;							//格子定数
		double kappa_s1, kappa_s2;						//勾配エネルギ－定数
		double ds_fac;										//phase fieldの揺らぎの大きさ

//---- 各種パラメータ設定 ----------------------------------------------------

		delt=0.1;	//時間きざみ入力

		temp=500.0;						//温度[K]
		al=250.0*1.0E-09;				//計算領域[m]
		b1=al/nd;							//差分ブロックの長さ[m]

		time1=-10.0;						//初期設定時間
		time1max=1.0+1.0e+07;		//計算時間の最大値

		smob=1.0;						//マルテンサイト変態ダイナミクスの緩和係数
		ds_fac=0.01;						//phase fieldの揺らぎ係数

		AA0=1000.0;						//マルテンサイト変態の化学的駆動力[J/mol]
		AA0=AA0/RR/temp;				//無次元化
		AA1=10.0; AA2=3.0*AA1+12.0; AA3=2.0*AA1+12.0;	//化学的駆動力定数

		kappa_s1=5.0e-15;												//勾配エネルギ－定数[Jm^2/mol]
		kappa_s1=kappa_s2=kappa_s1/RR/temp/b1/b1;			//無次元化

		a1_c=b1_c=c1_c=3.563e-10;									//格子定数[nm]
		atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//モル体積の計算（fccを仮定）[m^3/mol]

//--- s1場の変態歪の設定 ---
		eta_s1[1][1]=0.08; eta_s1[2][2]=-0.04; eta_s1[3][3]=0.0;
		eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.0;

//--- s2場の変態歪の設定 ---
		eta_s2[1][1]=eta_s1[2][2]; 	eta_s2[2][2]=eta_s1[1][1]; 	eta_s2[3][3]=0.0;
		eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.0;

//---弾性定数(Niの場合) ------------------------------------
		el_fac=1.0E+11*vm0/RR/temp;
		c11=2.508*el_fac;
		c44=1.235*el_fac;
		c12=1.500*el_fac;
		//c12=c11-2.0*c44;
		lam0=c12;	mu0=c44;			//ラーメの定数
		nu0=lam0/2.0/(lam0+mu0);	//ポアソン比

//---- 外力の設定 ---------------------------------------------
		sig22_a=-200.0*1.0e+06*vm0/RR/temp;						//ここでは外力を0に設定
		ep11_a=-lam0/4.0/mu0/(lam0+mu0)*sig22_a;					//平面歪を想定
		ep22_a=(lam0+2.0*mu0)/4.0/mu0/(lam0+mu0)*sig22_a;
		ep12_a=ep21_a=0.0;

//**** phase fieldの初期場設定と、sinおよびcosテーブルの設定 *******************************************
		prog.ini_field();				//phase fieldの初期場設定
		prog.table();					//フーリエ変換のためのsinとcosのテーブルとビット変換配列の設定

//**** phase fieldの時間発展の計算 *******************************************************************************
		while(time1<=time1max){

//---- phase fieldの表示 ----------------------------------------------------
			if((((int)(time1) % 50)==0)){ prog.update_draw(g); }		//カウント数が50の倍数おきに場を描画
			//if((((int)(time1) % 100)==0)){ prog.repaint(); } 

//---- phase fieldの保存 ----------------------------------------------------
			//if((((int)(time1) % 200)==0)){ prog.datsave(); }			//カウント数が200の倍数おきに場を保存
			//if(time1==3000.0){ prog.datsave(); }						//カウント数が3000の時に場を保存

//***** 勾配ポテンシャル *******************************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1; jp=j+1; jm=j-1;
					if(i==ndm){ip=0;}  if(i==0){im=ndm;}
					if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
					s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1h[i][j]);		//式(4)
					s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]-4.0*s2h[i][j]);		//式(4)
				}
			}

//**** 変態歪場のフ－リエ変換 ep11 ********************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]+eta_s2[1][1]*s2h[i][j];		//式(5)
					xi[i][j]=0.0;
				}
			}
			qs=-1.0; prog.rcfft();		//実空間からフーリエ空間への変換(qs<0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ 
					ep11qrh0[i][j]=xr[i][j];  ep11qih0[i][j]=xi[i][j];
				}
			}
			ep11qrh0[0][0]=ep11qih0[0][0]=0.0;

//**** 変態歪場のフ－リエ変換 ep22 ********************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]+eta_s2[2][2]*s2h[i][j]; 	//式(5)
					xi[i][j]=0.0;
				}
			}
			qs=-1.0; prog.rcfft();		//実空間からフーリエ空間への変換(qs<0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ 
					ep22qrh0[i][j]=xr[i][j];  ep22qih0[i][j]=xi[i][j];
				}
			}
			ep22qrh0[0][0]=ep22qih0[0][0]=0.0;

//*** 変態歪場の平均値の算出 ***
			sum11=sum22=0.0;
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j]; }
			}
			ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** 拘束歪変動量ec11の計算 *************************************
			for(i=0;i<=ndm;i++){
				if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
				for(j=0;j<=ndm;j++){
					if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
					alnn=Math.sqrt((double)ii*(double)ii+(double)jj*(double)jj);
					if(alnn==0.){alnn=1.;}
					nxx=(double)ii/alnn*(double)ii/alnn;
					nyy=(double)jj/alnn*(double)jj/alnn;
					Z11ep=nxx*(2.0*(1.0-nu0)-nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
					Z12ep=nxx*(2.0*nu0     -nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
					xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j]; 		//式(10)
					xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j]; 		//式(10)
				}
			}
			qs=1.0; prog.rcfft();		//フーリエ空間から実空間への変換(qs>0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ ec11[i][j]=xr[i][j]; }
			}

//***** 拘束歪変動量ec22の計算 *****************************
			for(i=0;i<=ndm;i++){
				if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
				for(j=0;j<=ndm;j++){
					if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
					alnn=Math.sqrt((double)ii*(double)ii+(double)jj*(double)jj);
					if(alnn==0.){alnn=1.;}
					nxx=(double)ii/alnn*(double)ii/alnn;
					nyy=(double)jj/alnn*(double)jj/alnn;
					Z21ep=nyy*(2.0*nu0     -nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
					Z22ep=nyy*(2.0*(1.0-nu0)-nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
					xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j]; 		//式(10)
					xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j]; 		//式(10)
				}
			}
			qs=1.0; prog.rcfft();		//フーリエ空間から実空間への変換(qs>0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ ec22[i][j]=xr[i][j]; }
			}

//****** ポテンシャルと発展方程式の計算 ************************************************************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					s1=s1h[i][j];  	s2=s2h[i][j];

//****** 化学ポテンシャルの計算 ********************************
					s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2)); 		//式(2)
					s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2)); 		//式(2)

//****** 弾性ポテンシャルの計算 ********************************
					ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
					ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;

					s1k_str=ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
								 +ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]); 		//式(7)
					s2k_str=ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
								 +ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]); 		//式(7)

//****** phase fieldの時間発展の計算 ********************************
					s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str); 				//式(12)
					s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str); 				//式(12)
					s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*Math.random()-1.0) )*delt;
					s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*Math.random()-1.0) )*delt;

//--- sの変域(0<=s<=1)の補正 ---
					if(s1h[i][j]>=1.0){s1h[i][j]=1.0;}  if(s1h[i][j]<=0.0){s1h[i][j]=0.0;}
					if(s2h[i][j]>=1.0){s2h[i][j]=1.0;}  if(s2h[i][j]<=0.0){s2h[i][j]=0.0;}
				}
			}

//**** [時間増加] *************************************************
			time1=time1+1.0;
		}//while

//----------------------------------------------------------------
		System.out.printf("\n 終了しました。グラフの右上×をクリックして終了してください。\n");
	}//main

// 以下はサブルーチンである。
// **** [phase fieldの初期設定] *****************************************************
	public void ini_field(){
		int i, j;
		double fac1;

		fac1=0.5; //最大の初期揺らぎ
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//s1h[i][j]=fac1*Math.random(); s2h[i][j]=fac1*Math.random();						//均一に核を置く場合
				if(Math.abs(j-nd2)<(nd/40)){s1h[i][j]=Math.random(); s2h[i][j]=Math.random();}	//中央に核を置く場合
			}
		}

	}

	//******* [フーリエ変換のためのsinとcosのテーブルと、ビット変換配列の設定] **********************
	public void table(){
		int it, it1, it2, mc, mn;
		double q;

		q=2.0*PI/(double)nd;
		for(it=0;it<=nd2-1;it++){ c[it]=Math.cos(q*(double)it); s[it]=Math.sin(q*(double)it); }

		ik[0]=0; mn=nd2; mc=1;
		for(it1=1;it1<=ig;it1++){
			for(it2=0;it2<=mc-1;it2++){ ik[it2+mc]=ik[it2]+mn; }
			mn=mn/2; mc=2*mc;
		}
	}

// *******************************************************************************
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** [phase fieldの描画] *******************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);//Windowをクリア

		int i, j, ii, jj;
		int icol, icol_r, icol_g, icol_b;
		double c_r, c_g, c_b;
		int ixmin=0, iymin=0, igx, igy, irad0;
		int ixmax=xwidth, iymax=yheight;
		double c, x, xmax, xmin, y, ymax, ymin, rad0;

		xmin=0.; xmax=1.; 						//横軸の最小値、最大値
		ymin=0.; ymax=1.; 						//縦軸の最小値、最大値
		rad0=1.0/(double)nd/2.0; 				//差分ブロックの長さの半分
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 ); 	//rad0のピクセル化

		System.out.printf("%f \n", time1); 	//計算の繰返し回数を標準入出力に表示

		for(i=0;i<=nd;i++){
			for(j=0;j<=nd;j++){
				//phase fieldの位置座標（実際の値）
				x=1.0/(double)nd*(double)i+rad0;
				y=1.0/(double)nd*(double)j+rad0;
				// phase fieldの位置座標（スクリーン座標に変換）
				igx=(int)( ((double)ixmax-(double)ixmin)*(x-xmin)/(xmax-xmin)+(double)ixmin );
				igy=(int)( ((double)iymax-(double)iymin)*(y-ymin)/(ymax-ymin)+(double)iymin );

				//個々の差分ブロックのphase field値
				ii=i; jj=j;
				if(i==nd){ii=0;} if(j==nd){jj=0;}			//周期的境界条件

				c_r= s1h[ii][jj];							//s1を赤
				c_g= s2h[ii][jj]; 							//s2を赤
				c_b=1.0-c_r-c_g; 						//変態前の相を青
				if(c_r>1.0){c_r=1.0;}  if(c_r<0.0){c_r=0.0;}
				if(c_g>1.0){c_g=1.0;}  if(c_g<0.0){c_g=0.0;}
				if(c_b>1.0){c_b=1.0;}  if(c_b<0.0){c_b=0.0;}

				icol_r=(int)(255.*c_r);  	icol_g=(int)(255.*c_g);  icol_b=(int)(255.*c_b);		//256階層に変換
				g.setColor(new Color(icol_r, icol_g, icol_b)); 									//差分ブロックの色を設定
				g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);					//個々の差分ブロック描画
			}
		}
	}

//***** [１次元高速フーリエ変換(FFT)] **************************************
	public void fft(){
		int ix, ka, kb, l2, lf, mf, n2, nf;
		double tj, tr;

		l2=1;
		for(lf=1;lf<=ig;lf++){
			n2=nd2/l2;
			for(mf=1;mf<=l2;mf++){
				for(nf=0;nf<=n2-1;nf++){
					ix=nf*l2; ka=nf+2*n2*(mf-1); kb=ka+n2;
					tr=xrf[ka]-xrf[kb];            tj=xif[ka]-xif[kb];
					xrf[ka]=xrf[ka]+xrf[kb];       xif[ka]=xif[ka]+xif[kb];
					xrf[kb]=tr*c[ix]+tj*qs*s[ix];    xif[kb]=tj*c[ix]-tr*qs*s[ix];
				}
			}
			l2=l2*2;
		}
	}

//**** [２次元高速フーリエ変換(RC FFT)] ***********************************
	public void rcfft(){
		int i, ic, ir, j;

		for(ir=0;ir<=ndm;ir++){
			for(ic=0;ic<=ndm;ic++){ xrf[ic]=xr[ir][ic];  xif[ic]=xi[ir][ic]; }
			fft();
			for(ic=0;ic<=ndm;ic++){ xr[ir][ic]=xrf[ik[ic]];  xi[ir][ic]=xif[ik[ic]]; }
		}
		for(ic=0;ic<=ndm;ic++){
			for(ir=0;ir<=ndm;ir++){ xrf[ir]=xr[ir][ic];  xif[ir]=xi[ir][ic]; }
			fft();
			for(ir=0;ir<=ndm;ir++){ xr[ir][ic]=xrf[ik[ir]];  xi[ir][ic]=xif[ik[ir]]; }
		}
		if(qs>0.0){return;}
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){ xr[i][j]=xr[i][j]/nd/nd;  xi[i][j]=xi[i][j]/nd/nd; }
		}
	}

//**** [デ－タの保存] ************************************
	private void datsave() throws Exception{
		int	i, j;

		//保存ファイル名をtest.datとする。
		PrintWriter outfile= new PrintWriter(							//ファイルのオープン
			new BufferedWriter(new FileWriter("test.dat", true)) );	//追記

		outfile.println(time1);											//カウントの書き込み
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				outfile.println(s1h[i][j]);									//phase fieldの書き込み
				outfile.println(s2h[i][j]);									//phase fieldの書き込み
			}
		}
		outfile.close();													//ファイルのクローズ
	}

//**** [デ－タの読込み] ************************************
	private void datin() throws Exception{
		int	i, j;
		String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//ファイルのオープン

		s_data=infile.readLine();  										//文字列として読み込み
		time1=new Double(s_data).doubleValue();					//文字を数値へ変換
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				s_data=infile.readLine();  								//文字列として読み込み
				s1h[i][j]=new Double(s_data).doubleValue();			//文字を数値へ変換
				s_data=infile.readLine(); 								//文字列として読み込み
				s2h[i][j]=new Double(s_data).doubleValue();			//文字を数値へ変換
			}
		}
		infile.close();														//ファイルのクローズ
	}

//*******************************************************************************************************
}//MT2D_001
//*** プログラム終了 ********************************************************************************
