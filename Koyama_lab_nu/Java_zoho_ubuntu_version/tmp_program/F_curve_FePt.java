//*** [プログラム (F_curve_FePt.java)] ************************************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class F_curve_FePt extends Frame{

	static int ND=201;		//濃度の分割数+1
	static int ndm=ND-1;	//濃度の分割数
	static int width;			// Window全体の幅
	static int height;			// Window全体の高さ
	static int xwidth;		// 描画領域の幅
	static int yheight;		// 描画領域の高さ
	static int insetx;			// Windowの枠の幅（左右および下）
	static int insety;			// Windowの枠の幅（上）
	static double PI=3.141592;					//π
	static double RR=8.3145;						//ガス定数
	static double [] sh=new double[ND];		//濃度
	static double [] Fh=new double[ND];		//自由エネルギー
	static double [] s2h=new double[ND];		//濃度
	static double [] F2h=new double[ND];		//自由エネルギー
	static int color_flg;							//(ch,Fh)を描くか、(c2h,F2h)を描くかを区別するフラグ
	static Graphics g;								//自由エネルギー曲線画面のグラフィックスオブジェクト

//*** [コンストラクタ] ****************************
	public F_curve_FePt(){
		xwidth=600; yheight=400; 		//描画画面の横と縦の長さ（ピクセル単位）
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

//*** [main プログラム] *****************************
	public static void main(String[] args) throws Exception{//例外処理は行わない

		F_curve_FePt prog=new F_curve_FePt();	//F_curve_FePtのインスタンスprogを生成

		int i;					//整数
		double s;		//規則度
		double AA1;	//係数

//---- AA1の値を設定 ----------------------------------------------------
		AA1=0.1;

//---- 標準入出力(コマンドプロンプト)より、係数値を入力する場合 --------
//    Javaの場合、標準入出力からの数値データ入力は、下記のように面倒であるので、
//    ここでは、参考までに定番の書式を掲載しておく。
//
//		String s_AA1;		//温度（文字列）
//		BufferedReader input=new BufferedReader(new InputStreamReader(System.in));
//		s_AA1="0.1";  //標準値を設定
//		try{ System.out.print("AA1(0.1) =  ");  s_AA1=input.readLine(); }
//		catch(IOException e){ System.out.println("Exception : "+e); }
//		AA1=new Double(s_AA1).doubleValue();		//文字列を数値に変換
//
//--- 自由エネルギーの計算 --------------------------------------------
		for(i=0;i<=ndm;i++){
			s=sh[i]=-1.0+2.0*(double)i/(double)ndm;
			Fh[i]=prog.G(s,AA1);  	//自由エネルギーのサブルーチンを呼んでいる。
			System.out.printf("%f   %e   \n", s, Fh[i]);//標準入出力へ自由エネルギー値を表示。
		}

//--- 自由エネルギーの描画 --------------------------------------
		color_flg=0; prog.repaint();

//--- データの保存 ---------------------------------------------
		prog.datsave();

//--- 描画が速すぎるので、強制的に5秒スリープ -----------------
		Thread.sleep(5000);
		//try{ Thread.sleep(5000); } catch( InterruptedException e){ }

//--- 上記で保存したデータを別の配列に読み込み -----------------
		prog.datin(); 

//--- 新しく読み込んだ配列の自由エネルギーを描画 --------------
		color_flg=1; prog.repaint();//再描画

//----------------------------------------------------------------

		System.out.printf("\n 終了しました。グラフの右上×をクリックして終了してください。\n");

}//main

// 以下はサブルーチンである。
// *** [自由エネルギー関数] ************************************************
	double G(double s, double AA1){
		double AA2, AA3;	//係数
		double gc;

		AA2=-4.0*AA1-12.0;  AA3=3.0*AA1+12.0;
		gc=0.5*AA1*s*s+0.25*AA2*s*s*s*s+1.0/6.0*AA3*s*s*s*s*s*s;

		return(gc);
	}

// *** [自由エネルギーのグラフ描画] **************************************
	public void paint(Graphics g){

		g.clearRect(0, 0, width, height);//Windowをクリア

		int i, i1, i2;
		double xmax, xmin, dx, ymax, ymin, dy;
		double gx1, gy1, gx2, gy2;
		int ixmax, iymax, ixmin, iymin, igx1, igy1, igx2, igy2;
		int idx, idy, ir;

		xmin=-1.0;     xmax=1.0;     dx=0.2;//横軸の最小値、最大値、分割間隔（実際の値）
		ymin=-1.0; ymax=0.2; dy=0.2;//縦軸の最小値、最大値、分割間隔（実際の値）
		ixmin=0; iymin=0;  ixmax=xwidth; iymax=yheight;//ピクセル単位の場合
		ir=4;//青丸の半径

		idx=(int)(ixmax*(dx/(xmax-xmin))+0.5);//ピクセル単位における横方向の分割間隔
		idy=(int)(iymax*(dy/(ymax-ymin))+0.5);//ピクセル単位における縦方向の分割間隔

		g.setColor(Color.white);  //色を白に指定
		g.fillRect(insetx, insety, xwidth, yheight);//画面を上で指定した色で塗る。

		g.setColor(Color.lightGray);//色を灰色に指定
		g.drawRect(insetx, insety, xwidth, yheight);//グラフの外側を上で指定した色で描く。
		for(i=0;i<=ixmax;i+=idx){ g.drawLine(insetx+i, insety+iymin, insetx+i, insety+iymax); }
		for(i=0;i<=iymax;i+=idy){ g.drawLine(insetx+ixmin, insety+i, insetx+ixmax, insety+i); }
				//さらにグラフ内の縦横の線を等間隔に描く。
//----------------------------------------------------------------------------------------------
		if(color_flg==0){//最初の描画
			g.setColor(Color.red);//色を赤に指定
			for(i=0;i<ndm;i++){
				i1=i; i2=i+1;			//個々の隣り合う２点を、直線で連続的に結ぶ。
				gx1=sh[i1];  gy1=Fh[i1];    gx2=sh[i2];  gy2=Fh[i2]; //実際の値

				igx1=(int)( ((double)ixmax-(double)ixmin)*(gx1-xmin)/(xmax-xmin)+(double)ixmin );
				igy1=(int)( (double)iymin+(double)iymax
						-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy1-ymin)+(double)iymin) );
				igx2=(int)( ((double)ixmax-(double)ixmin)*(gx2-xmin)/(xmax-xmin)+(double)ixmin );
				igy2=(int)( (double)iymin+(double)iymax
						-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy2-ymin)+(double)iymin) );
																								//スクリーン上のピクセル値に変換
   			g.drawLine(insetx+igx1, insety+igy1, insetx+igx2, insety+igy2);//２点を直線で結ぶ。
			}
		} 
		else{//再描画の場合
			g.setColor(Color.blue);//色を青に指定
			for(i=0;i<=ndm;i++){	//個々の点に小さな青丸を描く。
				gx1=s2h[i];  gy1=F2h[i]; //実際の値
				igx1=(int)( ((double)ixmax-(double)ixmin)*(gx1-xmin)/(xmax-xmin)+(double)ixmin );
				igy1=(int)( (double)iymin+(double)iymax
						-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy1-ymin)+(double)iymin) );
																								//スクリーン上のピクセル値に変換
				g.fillOval(insetx+igx1-ir, insety+igy1-ir, 2*ir, 2*ir);//青丸の描画
			}
		}

	}

//*** [デ－タの保存] ************************************
	private void datsave() throws Exception{
		int	i;

		//保存ファイル名をini000.datとする。
		PrintWriter outfile= new PrintWriter(
			new BufferedWriter(new FileWriter("ini000.dat")) );//ファイルオープン、上書きの場合
		//PrintWriter outfile= new PrintWriter(
			//new BufferedWriter(new FileWriter("ini000.dat", true)) );//ファイルオープン、追記の場合

		for(i=0;i<=ndm;i++){
			outfile.printf("%e  %e  \n", sh[i], Fh[i]);	//データの書き込み
			//outfile.println(sh[i]);	//データの書き込み
			//outfile.println(Fh[i]);	//データの書き込み
		}
		outfile.close();//ファイルのクローズ
	}

//*** [デ－タの読込み] ************************************
	private void datin() throws Exception{
		int	i;
		String s_data;
		String[] str_Ary;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//ファイルのオープン
			for(i=0;i<=ndm;i++){
				s_data=infile.readLine();  //文字列として読み込み
				str_Ary=s_data.split("  ");  //文字列を"  "にて２つに分割
				s2h[i]=new Double(str_Ary[0]).doubleValue();//文字を数値へ変換
				F2h[i]=new Double(str_Ary[1]).doubleValue();//文字を数値へ変換

				//s_data=infile.readLine();  //文字列をして読み込み
				//s2h[i]=new Double(s_data).doubleValue();//文字を数値へ変換
				//s_data=infile.readLine();   //文字列をして読み込み
				//F2h[i]=new Double(s_data).doubleValue();//文字を数値へ変換
			}
		infile.close();//ファイルのクローズ
	}

//****************************************************************
}//F_curve_FePt
//*** プログラム終了 ************************************************************

