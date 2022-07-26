import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class plot1D_c extends Frame{

	static int ND=512;		//�g�D�P�ӂ̕�����
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
	static double [] ch=new double[ND];	//�g�D���̔Z�x�f�|�^�z��
	static Graphics g;								//���R�G�l���M�[�Ȑ���ʂ̃O���t�B�b�N�X�I�u�W�F�N�g
	static double time1;	//�v�Z���ԁi�J�E���g���j

//*** [�R���X�g���N�^] ****************************
	public plot1D_c(){
		xwidth=800; yheight=200; 		//�`���ʂ̉��Əc�̒����i�s�N�Z���P�ʁj
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

		plot1D_c prog=new plot1D_c();

		int  i, j, k, l; 						//����
		double sumc;
    String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("c1D_field.dat"));
		while( (s_data=infile.readLine())!=null ){
			time1=new Double(s_data).doubleValue();
			for(i=0;i<=ndm;i++){
				s_data=infile.readLine();
				ch[i]=new Double(s_data).doubleValue();
			}

//--- �`�悪��������̂ŁA�����I��0.5�b�X���[�v -----------------
		Thread.sleep(500);

			prog.repaint();
		}
		infile.close();

		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");

	}//main

// **** �Z�x��`�� *****************************************************
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

		g.setColor(Color.black); 											//�`��F�����ɐݒ�
		g.setFont(new Font("Courier", Font.BOLD, 20));	//̫�Ă�ݒ�
		g.drawString("count="+time1, insetx, yheight+insety+25); 	//������`��

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