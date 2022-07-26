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
		jpg_file1="FePt_";  jpg_file0=jpg_file1;

//--- �ǂݏo���摜�t�@�C���̃J�E���g��(���R��) ------------------------
//--- �i���̒l���ۑ������摜�t�@�C�����ɒǉ������j ----------------
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
		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");

	}//main


// *********************************************************
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** [phase field�̕`��] *******************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);//Window���N���A

		int i, j, ii, jj;
		int icol, icol_r, icol_g, icol_b;
		double c1r, c1g, c1b, c2r, c2g, c2b;
		double c_r, c_g, c_b;
		int ixmin=0, iymin=0, igx, igy, irad0;
		int ixmax=xwidth, iymax=yheight;
		double c, x, xmax, xmin, y, ymax, ymin, rad0;

		xmin=0.; xmax=1.; 						//�����̍ŏ��l�A�ő�l
		ymin=0.; ymax=1.; 						//�c���̍ŏ��l�A�ő�l
		rad0=1.0/(double)nd/2.0; 				//�����u���b�N�̒����̔���
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 ); 	//rad0�̃s�N�Z����

		System.out.printf("%f \n", time1); 	//�v�Z�̌J�Ԃ��񐔂�W�����o�͂ɕ\��

		for(i=0;i<=nd;i++){
			for(j=0;j<=nd;j++){
				//phase field�̈ʒu���W�i���ۂ̒l�j
				x=1.0/(double)nd*(double)i+rad0;
				y=1.0/(double)nd*(double)j+rad0;
				// phase field�̈ʒu���W�i�X�N���[�����W�ɕϊ��j
				igx=(int)( ((double)ixmax-(double)ixmin)*(x-xmin)/(xmax-xmin)+(double)ixmin );
				igy=(int)( ((double)iymax-(double)iymin)*(y-ymin)/(ymax-ymin)+(double)iymin );

				//�X�̍����u���b�N��phase field�l
				ii=i; jj=j;
				if(i==nd){ii=0;} if(j==nd){jj=0;}			//�����I���E����

				c1r=c1g=c1b=c2r=c2g=c2b=0.0;
				if(s1h[ii][jj]>=0.0){c1r=s1h[ii][jj];c1g=0.0;c1b=0.0;}
				if(s1h[ii][jj]<0.0) {c1r=0.0;c1g=0.0;c1b=-s1h[ii][jj];}
				if(s2h[ii][jj]>=0.0){c2r=0.0;c2g=s2h[ii][jj];c2b=0.0;}
				if(s2h[ii][jj]<0.0) {c2r=-s2h[ii][jj];c2g=-s2h[ii][jj];c2b=0.0;}
				c_r=c1r+c2r;
				c_g=c1g+c2g;
				c_b=c1b+c2b;
				//c_r= s1h[ii][jj]*s1h[ii][jj];							//s1���
				//c_g= s2h[ii][jj]*s2h[ii][jj]; 						//s2���
				//c_b=1.0-c_r-c_g; 						//�ϑԑO�̑����

				if(c_r>1.0){c_r=1.0;}  if(c_r<0.0){c_r=0.0;}
				if(c_g>1.0){c_g=1.0;}  if(c_g<0.0){c_g=0.0;}
				if(c_b>1.0){c_b=1.0;}  if(c_b<0.0){c_b=0.0;}

				icol_r=(int)(255.*c_r);  	icol_g=(int)(255.*c_g);  icol_b=(int)(255.*c_b);	//256�K�w�ɕϊ�
				bg.setColor(new Color(icol_r, icol_g, icol_b)); 						//�����u���b�N�̐F��ݒ�
				bg.fillRect(igx-irad0,igy-irad0, irad0*2, irad0*2);	//�X�̍����u���b�N�`��
				//bg.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);	//�X�̍����u���b�N�`��

			}
		}
    g.drawImage(buff, 0, 0, this);
	}

//****************************************************************
}