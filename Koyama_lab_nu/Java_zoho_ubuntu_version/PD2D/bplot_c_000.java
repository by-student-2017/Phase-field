//表示される部分は上下左右が欠けるが、imageは全て保存される。

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.image.*;
import com.sun.image.codec.jpeg.*;

public class bplot_c extends Frame{

	static int ND=64;					//組織１辺の分割数
	static int nd=ND, ndm=ND-1;
	static int wh000, width, height;	//画面の幅と高さ
	static double [][] ch=new double[ND][ND];	//組織内の濃度デ－タ配列

	//グローバル変数化
	static double time1;	//計算時間（カウント数）
	//Image buff;
	static BufferedImage buff;
	static Graphics bg, g;

	//*************** コンストラクタ ****************************
	public bplot_c(){
		wh000=400; width=wh000;  height=wh000;  setSize(width, height);
		//wh000=400; width=wh000+4+4;  height=wh000+30+4;  setSize(width, height);
		setBackground(Color.lightGray);
		setVisible(true);
		//buff=this.createImage(width, height);
		buff=new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
		bg=buff.getGraphics();
		addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){ System.exit(0); }
		});
	}

//************** main *****************************
	public static void main(String[] args){

		bplot_c prog=new bplot_c();

		int  i, j, k, l; 						//整数
		double sumc;
    String str_data, str_time1;

		//prog.load("FeCu.jpg");

		try{
			BufferedReader infile=new BufferedReader(new FileReader("test.dat"));
			try{
				while( (str_time1=infile.readLine())!=null ){
					time1=new Double(str_time1).doubleValue();
					for(i=0;i<=ndm;i++){
						for(j=0;j<=ndm;j++){
							str_data=infile.readLine();
							ch[i][j]=new Double(str_data).doubleValue();
						}
					}
					prog.update_draw(g);
					if(time1==3000.0){prog.saveJpeg("FeCu.jpg");} //イメージをセーブするメソッド
				}
			}
			catch(Exception e){System.out.println(e); System.exit(1);}
			finally{infile.close();}
   	}
		catch(Exception e){System.out.println(e); System.exit(1);}

	}//main


// *********************************************************************
   public void update_draw( Graphics g ){ g=getGraphics(); draw(g); }

// **** 濃度場描画 *****************************************************
   public void draw(Graphics g){
      bg.clearRect(0, 0, width, height);
      bg.setColor(Color.lightGray);

			int xwidth=wh000, yheight=wh000, delxheight=0, delyheight=0;
			//int xwidth=wh000, yheight=wh000, delxheight=4, delyheight=30;
  	  int i, j, ii, jj;
    	int icol;
    	int ixmin=0, iymin=0, igx, igy, irad0;
			int ixmax=xwidth, iymax=yheight;
    	double c, x, xmax, xmin, y, ymax, ymin, rad0;

			xmin=0.; xmax=1.;  ymin=0.; ymax=1.;

			//printf("time %f\n",time1);
			System.out.println(time1);
			rad0=1.0/(double)nd/2.0;
			irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 );

			for(i=0;i<=nd;i++){
				for(j=0;j<=nd;j++){
					x=1.0/(double)nd*(double)i+rad0;
					igx=(int)( ((double)ixmax-(double)ixmin)*(x-xmin)/(xmax-xmin)+(double)ixmin );
					y=1.0/(double)nd*(double)j+rad0;
					igy=(int)( ((double)iymax-(double)iymin)*(y-ymin)/(ymax-ymin)+(double)iymin );
					ii=i; jj=j;
					if(i==nd){ii=0;}  if(j==nd){jj=0;}
					//icol=(int)(255.0*ch[ii][jj]);
					icol=(int)(255.0*(1.0-ch[ii][jj]));
					if(icol>=255){icol=255;}  if(icol<=0){icol=0;}
					bg.setColor(new Color(icol,icol,icol));
      		bg.fillRect(delxheight+igx-irad0,delyheight+igy-irad0, irad0*2, irad0*2);
					//bg.setColor(Color.blue);
					//bg.setColor(new Color(0,0,icol));
			}
		}
    g.drawImage(buff, 0, 0, this);
   }

//*** 画像の保存 **************************************************
		public void saveJpeg(String filename){
			try{
				FileOutputStream output=new FileOutputStream(filename);
				JPEGImageEncoder jpeg=JPEGCodec.createJPEGEncoder(output);
				jpeg.encode(buff);  output.flush();  output.close();
			}
			catch(IOException e){ e.printStackTrace(); }
		}

//*** 画像のロード **************************************************
		public void load(String fileName){
			try{
				FileInputStream input=new FileInputStream(fileName);
				JPEGImageDecoder decoder=JPEGCodec.createJPEGDecoder(input);
				buff=decoder.decodeAsBufferedImage();
				getGraphics().drawImage(buff,0,0,this);
			}
			catch(IOException e){}
		}

//****************************************************************
}