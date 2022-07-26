//表示される部分は上下左右が欠けるが、imageは全て保存される。

import javax.imageio.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.image.*;

public class bplot_c extends Frame{

	static int ND=64;		//組織１辺の分割数
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
	static double [][] ch=new double[ND][ND];	//組織内の濃度デ－タ配列
	static Graphics g, bg;								//自由エネルギー曲線画面のグラフィックスオブジェクト
	static double time1;	//計算時間（カウント数）
	static BufferedImage buff;

//*** [コンストラクタ] ****************************
	public bplot_c(){
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

		bplot_c prog=new bplot_c();

		int  i, j, k, l; 						//整数
		int  ii=1; 						//整数
		double sumc;
    String s_data;
    String jpg_file0, jpg_file1, jpg_file2;
		int [] jpg_n=new int[100];
		File outfile;

//--- 各画像ファイル名の共通部分 --------------------------------------
		jpg_file1="FeCr2D_";  jpg_file0=jpg_file1;

//--- 読み出す画像ファイルのカウント数(自然数) ------------------------
//--- （この値が保存される画像ファイル名に追加される） ----------------
		jpg_n[1]=0;
		jpg_n[2]=200;
		jpg_n[3]=600;
		jpg_n[4]=1000;
		jpg_n[5]=2000;
		jpg_n[6]=5000;
		jpg_n[7]=10000;
		jpg_n[8]=20000;
		//jpg_n[9]=3000;
		//jpg_n[10]=3000;
		//jpg_n[11]=3000;
		//jpg_n[12]=3000;


//----------------------------------------------------------------------
		BufferedReader infile=new BufferedReader(new FileReader("test_FeCr.dat"));
		while( (s_data=infile.readLine())!=null ){
			time1=new Double(s_data).doubleValue();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					s_data=infile.readLine();
					ch[i][j]=new Double(s_data).doubleValue();
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

	}//main


// *********************************************************
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** 濃度場描画 *****************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);//Windowをクリア

		int i, j, ii, jj;
		int icol;
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
				//icol=(int)(255.0*ch[ii][jj]);
				icol=(int)(255.0*(1.0-ch[ii][jj]));
				if(icol>=255){icol=255;} if(icol<=0){icol=0;}
				bg.setColor(new Color(icol,icol,icol));
      	bg.fillRect(igx-irad0,igy-irad0, irad0*2, irad0*2);//個々の差分ブロック描画
      	//bg.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);//個々の差分ブロック描画
			}
		}
    g.drawImage(buff, 0, 0, this);
	}

//****************************************************************
}