//*** [�v���O���� (PD1D_001.java)] ************************************************
//*** [�C���|�[�g��] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class PD1D_001 extends Frame{

//*** [�O���[�o���ϐ�] *****************************************************************************************
	static int ND=512;		//�g�D�P�ӂ̕�����
	static int nd=ND;			//�Z�x�̕�����
	static int ndm=ND-1;	//�Z�x�̕�����-1
	static int width;			// Window�S�̂̕�
	static int height;			// Window�S�̂̍���
	static int xwidth;			// �`��̈�̕�
	static int yheight;			// �`��̈�̍���
	static int insetx;			// Window�̘g�̕��i���E����щ��j
	static int insety;			// Window�̘g�̕��i��j
	static double PI=3.141592;				//��
	static double RR=8.3145;				//�K�X�萔
	static double [] ch=new double[ND];	//�g�D���̔Z�x�f�|�^�z��
	static Graphics g;						//���R�G�l���M�[�Ȑ���ʂ̃O���t�B�b�N�X�I�u�W�F�N�g
	static double time1;						//�v�Z���ԁi�J�E���g���j
	static double temp; 						//���x(K)
	static double c0;							//�����g���i���������j

//*** [�R���X�g���N�^] ****************************
	public PD1D_001(){
		xwidth=800; yheight=200; 			//�`���ʂ̉��Əc�̒����i�s�N�Z���P�ʁj
		insetx=4; insety=30;					//�`���ʂ̘g�̒���
		width=xwidth+insetx*2;  			//�`��Window�S�̂̉��̒���
		height=yheight+insetx+insety; 	//�`��Window�S�̂̏c�̒���
		setSize(width, height);				//�`��Window�̃Z�b�g
		setBackground(Color.white); 		//�`��Window�̕`�敔�̐F�𔒂ɐݒ�
		setVisible(true); 						//�`��Window��������悤�ɂ���
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
												//Window����鎞�̑���iWindow�̉E��~�̐ݒ�j
		});
	}

//*** [���C���v���O����] *******************************************************************************************
	public static void main(String[] args) throws Exception{//��O�����͍s��Ȃ�

		PD1D_001 prog=new PD1D_001();	//PD1D_001�̃C���X�^���Xprog�𐶐�

		int i;										//����
		int ip, im; 								//�����ii+1, i-1�j
		double delt;								//���ԍ��݁i�������j
		double al;									//�v�Z�̈�̂P�ӂ̒���
		double time1max;						//�v�Z���ԁi�ő�J�E���g���j
		double [] ck = new double[ND];		//�g�U�|�e���V����
		double [] ch2 = new double[ND];		//�g�D���̔Z�x�f�|�^�\���z��
		double mu_chem, mu_surf;			//�e�|�e���V����
		double c1, c2;							//�Z�x(A:1, B:2)
		double L0;								//���q�ԑ��ݍ�p�p�����[�^
		double kappa_c;							//�Z�x���z�G�l���M�[�W��
		double Mc;								//�Փ��x
		double c_flu;								//�Z�x��̗h�炬�̑傫��
		double cddtt;								//�Z�x�̑���

		double b1;								//�����u���b�N�T�C�Y
		double c2ip, c2im; 						//�����ɂ�����c2�𒆐S�ɁA���̍��E�̔Z�x
		double sumc, dc; 						//�Z�x��̑��a�A���ϑg������̂���

//---- �e��p�����[�^�ݒ� ----------------------------------------------------------------------
		temp=1000.0;			//[K]
		c0=0.1;					//�����̕��ϑg���i�����ł�A-40at%B������ݒ肵�Ă���j
		delt=0.01;
		time1=0.0;				//�v�Z�̌J�Ԃ���
		time1max=1.0e+08;		//�J�Ԃ��񐔂̍ő�l

		al=500.0;					//�Q�����v�Z�̈�̂P�ӂ̒���(nm)
		al=al*1.0e-9;				//(m)�ɕϊ�
		b1=al/(double)nd;		//�����P�u���b�N�̃T�C�Y

		L0=2.5e+04;				//���q�ԑ��ݍ�p�p�����[�^�iJ/mol�j
		L0=L0/RR/temp;			//��������

		kappa_c=5.0e-15;		//�Z�x���z�G�l���M�|�W���A�P�ʂ�[J m^2/mol]
		kappa_c=kappa_c/b1/b1/RR/temp;		//(b1^2*rr*temp)�Ŗ�������

		Mc=c0*(1.0-c0);			//�g�U�̈Փ��x
		c_flu=0.1;					//�Z�x��炬�U���̍ő�l

//---- ����0�ɂ����鏉���Z�x�v���t�@�C���ݒ� ------------------------------------------------------
		prog.ini_comp_field();	//�����ɂ���Đ�������ꍇ
		//prog.datin();				//�t�@�C������ǂݍ��ޏꍇ

//---- �Z�x�v���t�@�C���̎��Ԕ��W�̌v�Z ---------------------------------------------------------------
		while(time1<=time1max){

//---- �Z�x�v���t�@�C���̕\�� -------------------------------------------------------------------------------
			if((((int)(time1) % 500)==0)){ prog.repaint(); }	
											//�J�E���g����500�̔{�������ɔZ�x�v���t�@�C���`��

//---- �Z�x�v���t�@�C���̕ۑ� ------------------------------------------------------------------------------
			if((((int)(time1) % 1000)==0)){ prog.datsave(); }	
											//�J�E���g����1000�̔{�������ɔZ�x�v���t�@�C���ۑ�
			//if(time1==3000.0){ prog.datsave(); }	//�J�E���g����3000�̎��ɔZ�x�v���t�@�C���ۑ�

//---- �g�U�|�e���V�����̌v�Z -------------------------------------------------------------------------------
			for(i=0;i<=ndm;i++){
				ip=i+1; im=i-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}						//�����I���E����
				c2=ch[i];		c1=1.0-c2;									//�ʒui�ɂ�����c1��c2
				c2ip=ch[ip]; c2im=ch[im];									//�����ɂ�����c2�̍��E�̔Z�x

	 			mu_chem=L0*(c1-c2)+Math.log(c2)-Math.log(c1);		//���w�|�e���V������
				mu_surf=-2.0*kappa_c*(c2ip+c2im-2.0*c2);				//�Z�x���z�̃|�e���V����
				ck[i]=mu_chem+mu_surf; 									//�g�U�|�e���V����
			}

//---- �Z�x��̎��ԕω�(����`�g�U�������̍�����z��@) -------------------------------------------
			for(i=0;i<=ndm;i++){
				ip=i+1; im=i-1;  if(i==ndm){ip=0;}	if(i==0){im=ndm;} 	//�����I���E����
				cddtt=Mc*( ck[ip]+ ck[im]-2.0* ck[i]);						//����`�g�U������
				//ch2[i]=ch[i]+cddtt*delt; 									//�Z�x��̎��Ԕ��W
				ch2[i]=ch[i]+( cddtt+c_flu*(2.0*Math.random()-1.0) )*delt;	//�Z�x��̎��Ԕ��W�i�Z�x�h�炬���l���j
			}

//*** [�Z�x��̎��x�̕␳] *******************************************************
//*** ���l�v�Z�ł���̂ŔZ�x��̎��x�̕␳���s���i���ۂɂ͖��X�e�b�v�s���K�v�͂Ȃ��j****
  			sumc=0.; for(i=0;i<=ndm;i++){ sumc+=ch2[i]; }				//�Z�x�v���t�@�C���̐ϕ�
			dc=sumc/(double)nd-c0;										//�Z�x�v���t�@�C���̕ϓ���

			for(i=0;i<=ndm;i++){
				ch[i]=ch2[i]-dc;												//�Z�x��̕␳
				if(ch[i]>=1.){ch[i]=1.0-1.0e-6;}								//�Z�x��1�𒴂����ꍇ�̕␳
				if(ch[i]<=0.){ch[i]=1.0e-6;}									//�Z�x��0��؂����ꍇ�̕␳
			}

//******[���ԑ���]**************************************************************
			time1=time1+1.0;												//�v�Z���Ԃ̑���
		}//while

//----------------------------------------------------------------

		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");

}//main

// �ȉ��̓T�u���[�`��
// **** �����Z�x��ݒ� *****************************************************
	public void ini_comp_field(){
		int i;
		double rnd0, fac1;

		fac1=0.01; 				//�����Z�x�h�炬�̍ő�ω��ʂ�1%�ɐݒ�
		for(i=0;i<=ndm;i++){
			rnd0=2.0*Math.random()-1.0;  ch[i]=c0+rnd0*fac1;	//�����ɂď����Z�x���ݒ�
		}
	}

// **** �Z�x��`�� *****************************************************
	public void paint(Graphics g){
		int i, ii, i1, ii1, i2, ii2;
		double x, xmax, xmin, dx, y, ymax, ymin, dy, d0;
		double c, cmax, cmin;
		double gx1, gy1, gx2, gy2;
		int ixmin, iymin, ixmax, iymax, igx1, igy1, igx2, igy2, id0;
		int idx, idy;
		double col;

		xmin=0.; xmax=1.; dx=0.1; 								//�����̍ŏ��l�A�ő�l�A�����Ԋu�i���ۂ̒l�j
		ymin=0.; ymax=1.; dy=0.1;									//�c���̍ŏ��l�A�ő�l�A�����Ԋu�i���ۂ̒l�j
		cmin=0.0; cmax=1.0; 										//�Z�x�̍ŏ��l�A�ő�l�i���ۂ̒l�j
		ixmin=0; iymin=0;  ixmax=xwidth; iymax=yheight;		//�s�N�Z���P�ʂ̏ꍇ
		idx=(int)(0.1*ixmax);
		idy=(int)(0.1*iymax);

		g.setColor(Color.white); 									//�F�𔒂ɐݒ�
		g.fillRect(insetx, insety, xwidth, yheight); 				//��ʂ���Ŏw�肵���F�œh��
		g.setColor(Color.black); 									//�F�����ɐݒ�
		g.drawRect(insetx, insety, xwidth, yheight); 				//�O���t�̊O������Ŏw�肵���F�ŕ`��

		//�O���t���̏c���̐��𓙊Ԋu�ɕ`��
		for(i=0;i<=ixmax;i+=idx){ g.drawLine(insetx+i, insety+iymin, insetx+i, insety+iymax); }
		for(i=0;i<=iymax;i+=idy){ g.drawLine(insetx+ixmin, insety+i, insetx+ixmax, insety+i); }

		d0=1.0/(double)nd/2.0;										//�����u���b�N�̒����̔���
		id0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*d0 );		//d0�̃s�N�Z����

		System.out.printf("%f \n", time1);							//�v�Z�̌J�Ԃ��񐔂�W�����o�͂ɕ\��

		g.setColor(Color.red); 										//�F��Ԃɐݒ�
		for(i=0;i<=ndm;i++){
			//�Z�x�v���t�@�C���̒l�i���ۂ̒l�j
			i1=i; i2=i+1;
			gx1=1./(double)nd*(double)i1+d0;		gx2=1./(double)nd*(double)i2+d0;
			ii1=i1; ii2=i2;  if(i==ndm){ii2=0;}
			gy1=(ch[ii1]-cmin)/(cmax-cmin); 		gy2=(ch[ii2]-cmin)/(cmax-cmin);

			//�Z�x�v���t�@�C���̒l���X�N���[�����W�ɕϊ�
			igx1=(int)( ((double)ixmax-(double)ixmin)*(gx1-xmin)/(xmax-xmin)+(double)ixmin );
			igy1=(int)( (double)iymin+(double)iymax
									-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy1-ymin)+(double)iymin) );
			igx2=(int)( ((double)ixmax-(double)ixmin)*(gx2-xmin)/(xmax-xmin)+(double)ixmin );
			igy2=(int)( (double)iymin+(double)iymax
									-(((double)iymax-(double)iymin)/(ymax-ymin)*(gy2-ymin)+(double)iymin) );

			g.drawLine(insetx+igx1, insety+igy1, insetx+igx2, insety+igy2);	//�Z�x�v���t�@�C���̕`��
		}
	}

//*** [�f�|�^�̕ۑ�] ************************************
	private void datsave() throws Exception{
		int i;

		//�ۑ��t�@�C������test.dat�Ƃ���B
		PrintWriter outfile= new PrintWriter(
			new BufferedWriter(new FileWriter("test.dat", true)) );	//�t�@�C���̃I�[�v��(�ǋL)

		outfile.println(time1);								//�J�E���g�̏�������
		for(i=0;i<=ndm;i++){ outfile.println(ch[i]); }		//�Z�x��̏�������
		outfile.close();										//�t�@�C���̃N���[�Y
	}

//*** [�f�|�^�̓Ǎ���] ************************************
	private void datin() throws Exception{
		int i;
		String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//�t�@�C���̃I�[�v��

		s_data=infile.readLine();							//������Ƃ��ēǂݍ���
		time1=new Double(s_data).doubleValue();		//�����𐔒l�֕ϊ�
		for(i=0;i<=ndm;i++){
			s_data=infile.readLine();						//����������ēǂݍ���
			ch[i]=new Double(s_data).doubleValue();	//�����𐔒l�֕ϊ�
		}
		infile.close();											//�t�@�C���̃N���[�Y
	}

//****************************************************************
}//PD1D_001
//*** �v���O�����I�� ************************************************************
