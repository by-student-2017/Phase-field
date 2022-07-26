//�\������镔���͏㉺���E�������邪�Aimage�͑S�ĕۑ������B

import javax.imageio.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.image.*;

public class bplot_s12 extends Frame{

	static int ND=128;		//�g�D�P�ӂ̕�����
	static int nd=ND;		//�Z�x�̕�����
	static int ndm=ND-1;	//�Z�x�̕�����-1
	static int width;			// Window�S�̂̕�
	static int height;			// Window�S�̂̍���
	static int xwidth;		// �`��̈�̕�
	static int yheight;		// �`��̈�̍���
	static int insetx;			// Window�̘g�̕��i���E����щ��j
	static int insety;			// Window�̘g�̕��i��j
	static double PI=3.141592;					//��
	static double RR=8.3145;						//�K�X�萔
	static double [][] s1h=new double[ND][ND];	//�g�D���̔Z�x�f�|�^�z��
	static double [][] s2h=new double[ND][ND];	//�g�D���̔Z�x�f�|�^�z��
	static Graphics g, bg;								//���R�G�l���M�[�Ȑ���ʂ̃O���t�B�b�N�X�I�u�W�F�N�g
	static double time1;	//�v�Z���ԁi�J�E���g���j
	static BufferedImage buff;

//*** [�R���X�g���N�^] ****************************
	public bplot_s12(){
		xwidth=400; yheight=400; 		//�`���ʂ̉��Əc�̒����i�s�N�Z���P�ʁj
		insetx=4; insety=30;				//�`���ʂ̂ӂ��̒���
		width=xwidth;  							//
		height=yheight; 						//
		//width=xwidth+insetx*2;  		//�`��Window�S�̂̉��̒���
		//height=yheight+insetx+insety; 	//�`��Window�S�̂̏c�̒���
		setSize(width, height);			//�`��Window�̃Z�b�g
		setBackground(Color.white); 	//�`��Window�̕`�敔�̐F�𔒂ɐݒ�
		setVisible(true); 					//�`��Window��������悤�ɂ���
		buff=new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
		bg=buff.getGraphics();
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
												//Window����鎞�̑���iWindow�̉E��~�̐ݒ�j
		});
	}

//************** main *****************************
	public static void main(String[] args) throws Exception{

		bplot_s12 prog=new bplot_s12();

		int  i, j, k, l; 						//����
		int  ii=1; 						//����
		double sumc;
    String s_data;
    String jpg_file0, jpg_file1, jpg_file2;
		int [] jpg_n=new int[100];
		File outfile;

//--- �e�摜�t�@�C�����̋��ʕ��� --------------------------------------
		jpg_file1="MT2D_";  jpg_file0=jpg_file1;

//--- �ǂݏo���摜�t�@�C���̃J�E���g��(���R��) ------------------------
//--- �i���̒l���ۑ������摜�t�@�C�����ɒǉ������j ----------------
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
		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");

	}//main


// *********************************************************
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** ��`�� *****************************************************
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

				c_r= s1h[ii][jj];							//s1���
				c_g= s2h[ii][jj]; 							//s2���
				c_b=1.0-c_r-c_g; 						//�ϑԑO�̑����
				if(c_r>1.0){c_r=1.0;}  if(c_r<0.0){c_r=0.0;}
				if(c_g>1.0){c_g=1.0;}  if(c_g<0.0){c_g=0.0;}
				if(c_b>1.0){c_b=1.0;}  if(c_b<0.0){c_b=0.0;}

				icol_r=(int)(255.*c_r);  	icol_g=(int)(255.*c_g);  icol_b=(int)(255.*c_b);		//256�K�w�ɕϊ�
				bg.setColor(new Color(icol_r, icol_g, icol_b)); 									//�����u���b�N�̐F��ݒ�
      	bg.fillRect(igx-irad0,igy-irad0, irad0*2, irad0*2);//�X�̍����u���b�N�`��
				//bg.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);					//�X�̍����u���b�N�`��

			}
		}
    g.drawImage(buff, 0, 0, this);
	}

//****************************************************************
}