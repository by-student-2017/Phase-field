import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class plot_c extends Frame{

	static int ND=64;		//�g�D�P�ӂ̕�����
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
	static double [][] ch=new double[ND][ND];	//�g�D���̔Z�x�f�|�^�z��
	static Graphics g;								//���R�G�l���M�[�Ȑ���ʂ̃O���t�B�b�N�X�I�u�W�F�N�g
	static double time1;	//�v�Z���ԁi�J�E���g���j

//*** [�R���X�g���N�^] ****************************
	public plot_c(){
		xwidth=400; yheight=400; 		//�`���ʂ̉��Əc�̒����i�s�N�Z���P�ʁj
		insetx=4; insety=30;				//�`���ʂ̂ӂ��̒���
		width=xwidth+insetx*2;  		//�`��Window�S�̂̉��̒���
		height=yheight+insetx+insety; 	//�`��Window�S�̂̏c�̒���
		setSize(width, height);			//�`��Window�̃Z�b�g
		setBackground(Color.white); 	//�`��Window�̕`�敔�̐F�𔒂ɐݒ�
		setVisible(true); 					//�`��Window��������悤�ɂ���
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
												//Window����鎞�̑���iWindow�̉E��~�̐ݒ�j
		});
	}

//************** main *****************************
	public static void main(String[] args) throws Exception{

		plot_c prog=new plot_c();

		int  i, j, k, l; 						//����
		double sumc;
    String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("test_FeCr.dat"));
		while( (s_data=infile.readLine())!=null ){
			time1=new Double(s_data).doubleValue();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					s_data=infile.readLine();
					ch[i][j]=new Double(s_data).doubleValue();
				}
			}
			prog.update_draw(g);
		}
		infile.close();

		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");

	}//main

// *********************************************************
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** �Z�x��`�� *****************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);//Window���N���A

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
				g.setColor(new Color(icol,icol,icol));
      	g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);//�X�̍����u���b�N�`��
			}
		}
	}

//****************************************************************
}