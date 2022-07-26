//表示される部分は上下左右が欠けるが、imageは全て保存される。

import javax.imageio.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.image.*;

public class bplot_s12 extends Frame{

	static int ND=128;		//組織１辺の分割数
	static int nd=ND;		//濃度の分割数
	static int ndm=ND-1;	//濃度の分割数-1
	static int width;			// Window全体の幅
	static int height;			// Window全体の高さ
	static int xwidth;		// 描画領域の幅
	static int yheight;		// 描画領域の高さ
	static int insetx;			// Windowの枠の幅（左右および下）
	static int insety;			// Windowの枠の幅（上）
	static double PI=3.141592;					//π
	static double RR=8.3145;						//ガス定数
	static double [][] s1h=new double[ND][ND];	//組織内の濃度デ－タ配列
	static double [][] s2h=new double[ND][ND];	//組織内の濃度デ－タ配列
	static Graphics g, bg;								//自由エネルギー曲線画面のグラフィックスオブジェクト
	static double time1;	//計算時間（カウント数）
	static BufferedImage buff;

//*** [コンストラクタ] ****************************
	public bplot_s12(){
		xwidth=400; yheight=400; 		//描画画面の横と縦の長さ（ピクセル単位）
		insetx=4; insety=30;				//描画画面のふちの長さ
		width=xwidth;  							//
		height=yheight; 						//
		//width=xwidth+insetx*2;  		//描画Window全体の横の長さ
		//height=yheight+insetx+insety; 	//描画Window全体の縦の長さ
		setSize(width, height);			//描画Windowのセット
		setBackground(Color.white); 	//描画Windowの描画部の色を白に設定
		setVisible(true); 					//描画Windowを見えるようにする
		buff=new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
		bg=buff.getGraphics();
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
												//Windowを閉じる時の操作（Windowの右上×の設定）
		});
	}

//************** main *****************************
	public static void main(String[] args) throws Exception{

		bplot_s12 prog=new bplot_s12();

		int  i, j, k, l; 						//整数
		int  ii=1; 						//整数
		double sumc;
    String s_data;
    String jpg_file0, jpg_file1, jpg_file2;
		int [] jpg_n=new int[100];
		File outfile;

//--- 各画像ファイル名の共通部分 --------------------------------------
		jpg_file1="FePt_";  jpg_file0=jpg_file1;

//--- 読み出す画像ファイルのカウント数(自然数) ------------------------
//--- （この値が保存される画像ファイル名に追加される） ----------------
		jpg_n[1]=0;
		jpg_n[2]=100;
		jpg_n[3]=200;
		jpg_n[4]=300;
		jpg_n[5]=400;
		jpg_n[6]=500;
		jpg_n[7]=600;
		jpg_n[8]=700;
		jpg_n[9]=800;
		jpg_n[10]=900;
		jpg_n[11]=1000;
		jpg_n[12]=2000;
		jpg_n[13]=3000;
		jpg_n[14]=4000;
		jpg_n[15]=5000;
		jpg_n[16]=10000;
		//jpg_n[15]=3000;
		//jpg_n[16]=3000;


//----------------------------------------------------------------------
		BufferedReader infile=new BufferedReader(new FileReader("test_FePt.dat"));
		while( (s_data=infile.readLine())!=null ){
			time1=new Double(s_data).doubleValue();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					s_data=infile.readLine();
					s1h[i][j]=new Double(s_data).doubleValue();
					s_data=infile.readLine();
					s2h[i][j]=new Double(s_data).doubleValue();
				}
			}
			if((int)time1==jpg_n[ii]){
				prog.update_draw(g);
				jpg_file2=Integer.toString(jpg_n[ii]);  
				jpg_file1=jpg_file1+jpg_file2+".jpg";
				outfile=new File(jpg_file1);  
				ImageIO.write(buff, "jpeg", outfile);
				ii=ii+1;  jpg_file1=jpg_file0;
			}
		}
		infile.close();

//----------------------------------------------------------------
		System.out.printf("\n 終了しました。グラフの右上×をクリックして終了してください。\n");

	}//main


// *********************************************************
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** [phase fieldの描画] *******************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);//Windowをクリア

		int i, j, ii, jj;
		int icol, icol_r, icol_g, icol_b;
		double c1r, c1g, c1b, c2r, c2g, c2b;
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

				c1r=c1g=c1b=c2r=c2g=c2b=0.0;
				if(s1h[ii][jj]>=0.0){c1r=s1h[ii][jj];c1g=0.0;c1b=0.0;}
				if(s1h[ii][jj]<0.0) {c1r=0.0;c1g=0.0;c1b=-s1h[ii][jj];}
				if(s2h[ii][jj]>=0.0){c2r=0.0;c2g=s2h[ii][jj];c2b=0.0;}
				if(s2h[ii][jj]<0.0) {c2r=-s2h[ii][jj];c2g=-s2h[ii][jj];c2b=0.0;}
				c_r=c1r+c2r;
				c_g=c1g+c2g;
				c_b=c1b+c2b;
				//c_r= s1h[ii][jj]*s1h[ii][jj];							//s1を赤
				//c_g= s2h[ii][jj]*s2h[ii][jj]; 						//s2を緑
				//c_b=1.0-c_r-c_g; 						//変態前の相を青

				if(c_r>1.0){c_r=1.0;}  if(c_r<0.0){c_r=0.0;}
				if(c_g>1.0){c_g=1.0;}  if(c_g<0.0){c_g=0.0;}
				if(c_b>1.0){c_b=1.0;}  if(c_b<0.0){c_b=0.0;}

				icol_r=(int)(255.*c_r);  	icol_g=(int)(255.*c_g);  icol_b=(int)(255.*c_b);	//256階層に変換
				bg.setColor(new Color(icol_r, icol_g, icol_b)); 						//差分ブロックの色を設定
				bg.fillRect(igx-irad0,igy-irad0, irad0*2, irad0*2);	//個々の差分ブロック描画
				//bg.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);	//個々の差分ブロック描画

			}
		}
    g.drawImage(buff, 0, 0, this);
	}

//****************************************************************
}