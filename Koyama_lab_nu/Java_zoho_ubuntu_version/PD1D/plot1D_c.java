import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class plot1D_c extends Frame{

	static int ND=512;		//組織１辺の分割数
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
	static double [] ch=new double[ND];	//組織内の濃度デ−タ配列
	static Graphics g;								//自由エネルギー曲線画面のグラフィックスオブジェクト
	static double time1;	//計算時間（カウント数）

//*** [コンストラクタ] ****************************
	public plot1D_c(){
		xwidth=800; yheight=200; 		//描画画面の横と縦の長さ（ピクセル単位）
		insetx=4; insety=30;				//描画画面のふちの長さ
		width=xwidth+insetx*2;  		//描画Window全体の横の長さ
		height=yheight+insetx+insety; 	//描画Window全体の縦の長さ
		setSize(width, height);			//描画Windowのセット
		setBackground(Color.white); 	//描画Windowの描画部の色を白に設定
		setVisible(true); 					//描画Windowを見えるようにする
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
												//Windowを閉じる時の操作（Windowの右上×の設定）
		});
	}

//************** main *****************************
	public static void main(String[] args) throws Exception{

		plot1D_c prog=new plot1D_c();

		int  i, j, k, l; 						//整数
		double sumc;
    String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("c1D_field.dat"));
		while( (s_data=infile.readLine())!=null ){
			time1=new Double(s_data).doubleValue();
			for(i=0;i<=ndm;i++){
				s_data=infile.readLine();
				ch[i]=new Double(s_data).doubleValue();
			}

//--- 描画が速すぎるので、強制的に0.5秒スリープ -----------------
		Thread.sleep(500);

			prog.repaint();
		}
		infile.close();

		System.out.printf("\n 終了しました。グラフの右上×をクリックして終了してください。\n");

	}//main

// **** 濃度場描画 *****************************************************
	public void paint(Graphics g){
		int ixmax=xwidth, iymax=yheight;
		int i, ii, i1, ii1, i2, ii2;
		double col;
		int ixmin=0, iymin=0, igx1, igy1, igx2, igy2, irad0;
		double cmax, cmin, dx, dy;
		int idx, idy;
    double c, x, xmax, xmin, y, ymax, ymin, rad0;
		double gx1, gy1, gx2, gy2;

		xmin=0.; xmax=1.;  ymin=0.; ymax=1.;
		cmax=1.0; cmin=0.0;
		dx=0.1; dy=0.1;
		idx=(int)(0.1*ixmax);
		idy=(int)(0.1*iymax);

		g.setColor(Color.black); 											//描画色を黒に設定
		g.setFont(new Font("Courier", Font.BOLD, 20));	//ﾌｫﾝﾄを設定
		g.drawString("count="+time1, insetx, yheight+insety+25); 	//文字列描画

		g.setColor(Color.white);
   	g.fillRect(insetx, insety, xwidth, yheight);
		g.setColor(Color.black);
   	g.drawRect(insetx, insety, xwidth, yheight);

		for(i=0;i<=ixmax;i+=idx){
    	g.drawLine(insetx+i, insety+iymin, insetx+i, insety+iymax);
		}
		for(i=0;i<=iymax;i+=idy){
    	g.drawLine(insetx+ixmin, insety+i, insetx+ixmax, insety+i);
		}

		rad0=1.0/(double)nd/2.0;
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 );

		System.out.printf("%f \n", time1);

		g.setColor(Color.red);
		for(i=0;i<=ndm;i++){
			i1=i; i2=i+1;
			gx1=1./(double)nd*(double)i1+rad0;
			gx2=1./(double)nd*(double)i2+rad0;
			ii1=i1; ii2=i2;  if(i==ndm){ii2=0;}

			gy1=(ch[ii1]-cmin)/(cmax-cmin);
			gy2=(ch[ii2]-cmin)/(cmax-cmin);

			igx1=(int)( ((double)ixmax-(double)ixmin)*(gx1-xmin)/(xmax-xmin)+(double)ixmin );
			igy1=(int)( (double)iymin+(double)iymax
									-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy1-ymin)+(double)iymin) );
			igx2=(int)( ((double)ixmax-(double)ixmin)*(gx2-xmin)/(xmax-xmin)+(double)ixmin );
			igy2=(int)( (double)iymin+(double)iymax
									-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy2-ymin)+(double)iymin) );

    	g.drawLine(insetx+igx1, insety+igy1, insetx+igx2, insety+igy2);
		}
	}

//****************************************************************
}