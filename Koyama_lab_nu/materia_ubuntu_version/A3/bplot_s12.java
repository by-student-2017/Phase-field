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
		jpg_file1="MT2D_";  jpg_file0=jpg_file1;

//--- 読み出す画像ファイルのカウント数(自然数) ------------------------
//--- （この値が保存される画像ファイル名に追加される） ----------------
		jpg_n[1]=0;
		jpg_n[2]=200;
		jpg_n[3]=400;
		jpg_n[4]=600;
		jpg_n[5]=800;
		jpg_n[6]=1000;
		jpg_n[7]=1200;
		jpg_n[8]=1400;
		jpg_n[9]=1600;
		jpg_n[10]=1800;
		jpg_n[11]=2000;
		//jpg_n[12]=3000;


//----------------------------------------------------------------------
		BufferedReader infile=new BufferedReader(new FileReader("test.dat"));
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

// **** 場描画 *****************************************************
	public void paint(Graphics g){

		int i, j, ii, jj;
		int icol, icol_r, icol_g, icol_b;
		double c_r, c_g, c_b;
		int ixmin=0, iymin=0, igx, igy, irad0;
		int ixmax=xwidth, iymax=yheight;
		double c, x, xmax, xmin, y, ymax, ymin, rad0;

		xmin=0.; xmax=1.;  ymin=0.; ymax=1.;

		rad0=1.0/(double)nd/2.0;
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 );

		System.out.printf("%f \n", time1);

		for(i=0;i<=nd;i++){
			for(j=0;j<=nd;j++){
				x=1.0/(double)nd*(double)i+rad0;
				igx=(int)( ((double)ixmax-(double)ixmin)*(x-xmin)/(xmax-xmin)+(double)ixmin );
				y=1.0/(double)nd*(double)j+rad0;
				igy=(int)( ((double)iymax-(double)iymin)*(y-ymin)/(ymax-ymin)+(double)iymin );
				ii=i; jj=j;
				if(i==nd){ii=0;} if(j==nd){jj=0;}

				c_r= s1h[ii][jj];							//s1を赤
				c_g= s2h[ii][jj]; 							//s2を赤
				c_b=1.0-c_r-c_g; 						//変態前の相を青
				if(c_r>1.0){c_r=1.0;}  if(c_r<0.0){c_r=0.0;}
				if(c_g>1.0){c_g=1.0;}  if(c_g<0.0){c_g=0.0;}
				if(c_b>1.0){c_b=1.0;}  if(c_b<0.0){c_b=0.0;}

				icol_r=(int)(255.*c_r);  	icol_g=(int)(255.*c_g);  icol_b=(int)(255.*c_b);		//256階層に変換
				bg.setColor(new Color(icol_r, icol_g, icol_b)); 									//差分ブロックの色を設定
      	bg.fillRect(igx-irad0,igy-irad0, irad0*2, irad0*2);//個々の差分ブロック描画
				//bg.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);					//個々の差分ブロック描画

			}
		}
    g.drawImage(buff, 0, 0, this);
	}

//****************************************************************
}