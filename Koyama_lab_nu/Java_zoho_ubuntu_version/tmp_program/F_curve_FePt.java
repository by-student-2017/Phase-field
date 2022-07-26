//*** [�v���O���� (F_curve_FePt.java)] ************************************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class F_curve_FePt extends Frame{

	static int ND=201;		//�Z�x�̕�����+1
	static int ndm=ND-1;	//�Z�x�̕�����
	static int width;			// Window�S�̂̕�
	static int height;			// Window�S�̂̍���
	static int xwidth;		// �`��̈�̕�
	static int yheight;		// �`��̈�̍���
	static int insetx;			// Window�̘g�̕��i���E����щ��j
	static int insety;			// Window�̘g�̕��i��j
	static double PI=3.141592;					//��
	static double RR=8.3145;						//�K�X�萔
	static double [] sh=new double[ND];		//�Z�x
	static double [] Fh=new double[ND];		//���R�G�l���M�[
	static double [] s2h=new double[ND];		//�Z�x
	static double [] F2h=new double[ND];		//���R�G�l���M�[
	static int color_flg;							//(ch,Fh)��`�����A(c2h,F2h)��`��������ʂ���t���O
	static Graphics g;								//���R�G�l���M�[�Ȑ���ʂ̃O���t�B�b�N�X�I�u�W�F�N�g

//*** [�R���X�g���N�^] ****************************
	public F_curve_FePt(){
		xwidth=600; yheight=400; 		//�`���ʂ̉��Əc�̒����i�s�N�Z���P�ʁj
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

//*** [main �v���O����] *****************************
	public static void main(String[] args) throws Exception{//��O�����͍s��Ȃ�

		F_curve_FePt prog=new F_curve_FePt();	//F_curve_FePt�̃C���X�^���Xprog�𐶐�

		int i;					//����
		double s;		//�K���x
		double AA1;	//�W��

//---- AA1�̒l��ݒ� ----------------------------------------------------
		AA1=0.1;

//---- �W�����o��(�R�}���h�v�����v�g)���A�W���l����͂���ꍇ --------
//    Java�̏ꍇ�A�W�����o�͂���̐��l�f�[�^���͂́A���L�̂悤�ɖʓ|�ł���̂ŁA
//    �����ł́A�Q�l�܂łɒ�Ԃ̏������f�ڂ��Ă����B
//
//		String s_AA1;		//���x�i������j
//		BufferedReader input=new BufferedReader(new InputStreamReader(System.in));
//		s_AA1="0.1";  //�W���l��ݒ�
//		try{ System.out.print("AA1(0.1) =  ");  s_AA1=input.readLine(); }
//		catch(IOException e){ System.out.println("Exception : "+e); }
//		AA1=new Double(s_AA1).doubleValue();		//������𐔒l�ɕϊ�
//
//--- ���R�G�l���M�[�̌v�Z --------------------------------------------
		for(i=0;i<=ndm;i++){
			s=sh[i]=-1.0+2.0*(double)i/(double)ndm;
			Fh[i]=prog.G(s,AA1);  	//���R�G�l���M�[�̃T�u���[�`�����Ă�ł���B
			System.out.printf("%f   %e   \n", s, Fh[i]);//�W�����o�͂֎��R�G�l���M�[�l��\���B
		}

//--- ���R�G�l���M�[�̕`�� --------------------------------------
		color_flg=0; prog.repaint();

//--- �f�[�^�̕ۑ� ---------------------------------------------
		prog.datsave();

//--- �`�悪��������̂ŁA�����I��5�b�X���[�v -----------------
		Thread.sleep(5000);
		//try{ Thread.sleep(5000); } catch( InterruptedException e){ }

//--- ��L�ŕۑ������f�[�^��ʂ̔z��ɓǂݍ��� -----------------
		prog.datin(); 

//--- �V�����ǂݍ��񂾔z��̎��R�G�l���M�[��`�� --------------
		color_flg=1; prog.repaint();//�ĕ`��

//----------------------------------------------------------------

		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");

}//main

// �ȉ��̓T�u���[�`���ł���B
// *** [���R�G�l���M�[�֐�] ************************************************
	double G(double s, double AA1){
		double AA2, AA3;	//�W��
		double gc;

		AA2=-4.0*AA1-12.0;  AA3=3.0*AA1+12.0;
		gc=0.5*AA1*s*s+0.25*AA2*s*s*s*s+1.0/6.0*AA3*s*s*s*s*s*s;

		return(gc);
	}

// *** [���R�G�l���M�[�̃O���t�`��] **************************************
	public void paint(Graphics g){

		g.clearRect(0, 0, width, height);//Window���N���A

		int i, i1, i2;
		double xmax, xmin, dx, ymax, ymin, dy;
		double gx1, gy1, gx2, gy2;
		int ixmax, iymax, ixmin, iymin, igx1, igy1, igx2, igy2;
		int idx, idy, ir;

		xmin=-1.0;     xmax=1.0;     dx=0.2;//�����̍ŏ��l�A�ő�l�A�����Ԋu�i���ۂ̒l�j
		ymin=-1.0; ymax=0.2; dy=0.2;//�c���̍ŏ��l�A�ő�l�A�����Ԋu�i���ۂ̒l�j
		ixmin=0; iymin=0;  ixmax=xwidth; iymax=yheight;//�s�N�Z���P�ʂ̏ꍇ
		ir=4;//�ۂ̔��a

		idx=(int)(ixmax*(dx/(xmax-xmin))+0.5);//�s�N�Z���P�ʂɂ����鉡�����̕����Ԋu
		idy=(int)(iymax*(dy/(ymax-ymin))+0.5);//�s�N�Z���P�ʂɂ�����c�����̕����Ԋu

		g.setColor(Color.white);  //�F�𔒂Ɏw��
		g.fillRect(insetx, insety, xwidth, yheight);//��ʂ���Ŏw�肵���F�œh��B

		g.setColor(Color.lightGray);//�F���D�F�Ɏw��
		g.drawRect(insetx, insety, xwidth, yheight);//�O���t�̊O������Ŏw�肵���F�ŕ`���B
		for(i=0;i<=ixmax;i+=idx){ g.drawLine(insetx+i, insety+iymin, insetx+i, insety+iymax); }
		for(i=0;i<=iymax;i+=idy){ g.drawLine(insetx+ixmin, insety+i, insetx+ixmax, insety+i); }
				//����ɃO���t���̏c���̐��𓙊Ԋu�ɕ`���B
//----------------------------------------------------------------------------------------------
		if(color_flg==0){//�ŏ��̕`��
			g.setColor(Color.red);//�F��ԂɎw��
			for(i=0;i<ndm;i++){
				i1=i; i2=i+1;			//�X�ׂ̗荇���Q�_���A�����ŘA���I�Ɍ��ԁB
				gx1=sh[i1];  gy1=Fh[i1];    gx2=sh[i2];  gy2=Fh[i2]; //���ۂ̒l

				igx1=(int)( ((double)ixmax-(double)ixmin)*(gx1-xmin)/(xmax-xmin)+(double)ixmin );
				igy1=(int)( (double)iymin+(double)iymax
						-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy1-ymin)+(double)iymin) );
				igx2=(int)( ((double)ixmax-(double)ixmin)*(gx2-xmin)/(xmax-xmin)+(double)ixmin );
				igy2=(int)( (double)iymin+(double)iymax
						-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy2-ymin)+(double)iymin) );
																								//�X�N���[����̃s�N�Z���l�ɕϊ�
   			g.drawLine(insetx+igx1, insety+igy1, insetx+igx2, insety+igy2);//�Q�_�𒼐��Ō��ԁB
			}
		} 
		else{//�ĕ`��̏ꍇ
			g.setColor(Color.blue);//�F��Ɏw��
			for(i=0;i<=ndm;i++){	//�X�̓_�ɏ����Ȑۂ�`���B
				gx1=s2h[i];  gy1=F2h[i]; //���ۂ̒l
				igx1=(int)( ((double)ixmax-(double)ixmin)*(gx1-xmin)/(xmax-xmin)+(double)ixmin );
				igy1=(int)( (double)iymin+(double)iymax
						-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy1-ymin)+(double)iymin) );
																								//�X�N���[����̃s�N�Z���l�ɕϊ�
				g.fillOval(insetx+igx1-ir, insety+igy1-ir, 2*ir, 2*ir);//�ۂ̕`��
			}
		}

	}

//*** [�f�|�^�̕ۑ�] ************************************
	private void datsave() throws Exception{
		int	i;

		//�ۑ��t�@�C������ini000.dat�Ƃ���B
		PrintWriter outfile= new PrintWriter(
			new BufferedWriter(new FileWriter("ini000.dat")) );//�t�@�C���I�[�v���A�㏑���̏ꍇ
		//PrintWriter outfile= new PrintWriter(
			//new BufferedWriter(new FileWriter("ini000.dat", true)) );//�t�@�C���I�[�v���A�ǋL�̏ꍇ

		for(i=0;i<=ndm;i++){
			outfile.printf("%e  %e  \n", sh[i], Fh[i]);	//�f�[�^�̏�������
			//outfile.println(sh[i]);	//�f�[�^�̏�������
			//outfile.println(Fh[i]);	//�f�[�^�̏�������
		}
		outfile.close();//�t�@�C���̃N���[�Y
	}

//*** [�f�|�^�̓Ǎ���] ************************************
	private void datin() throws Exception{
		int	i;
		String s_data;
		String[] str_Ary;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//�t�@�C���̃I�[�v��
			for(i=0;i<=ndm;i++){
				s_data=infile.readLine();  //������Ƃ��ēǂݍ���
				str_Ary=s_data.split("  ");  //�������"  "�ɂĂQ�ɕ���
				s2h[i]=new Double(str_Ary[0]).doubleValue();//�����𐔒l�֕ϊ�
				F2h[i]=new Double(str_Ary[1]).doubleValue();//�����𐔒l�֕ϊ�

				//s_data=infile.readLine();  //����������ēǂݍ���
				//s2h[i]=new Double(s_data).doubleValue();//�����𐔒l�֕ϊ�
				//s_data=infile.readLine();   //����������ēǂݍ���
				//F2h[i]=new Double(s_data).doubleValue();//�����𐔒l�֕ϊ�
			}
		infile.close();//�t�@�C���̃N���[�Y
	}

//****************************************************************
}//F_curve_FePt
//*** �v���O�����I�� ************************************************************

