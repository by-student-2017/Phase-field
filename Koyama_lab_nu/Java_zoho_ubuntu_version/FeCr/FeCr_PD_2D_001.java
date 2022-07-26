//*** [�v���O���� (FeCr_PD_2D_001.java)] ************************************************
//*** [�C���|�[�g��] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class FeCr_PD_2D_001 extends Frame{

//*** [�O���[�o���ϐ�] *****************************************************************************************
	static int ND=64;			//�g�D�P�ӂ̕�����
	static int nd=ND;			//�Z�x�̕�����
	static int ndm=ND-1;	//�Z�x�̕�����-1
	static int width;			// Window�S�̂̕�
	static int height;			// Window�S�̂̍���
	static int xwidth;			// �`��̈�̕�
	static int yheight;			// �`��̈�̍���
	static int insetx;			// Window�̘g�̕��i���E����щ��j
	static int insety;			// Window�̘g�̕��i��j
	static double PI=3.141592;						//��
	static double RR=8.3145;						//�K�X�萔
	static double [][] ch=new double[ND][ND];	//�g�D���̔Z�x�f�|�^�z��
	static Graphics g;								//���R�G�l���M�[�Ȑ���ʂ̃O���t�B�b�N�X�I�u�W�F�N�g
	static double time1;								//�v�Z���ԁi�J�E���g���j
	static double temp; 								//���x(K)
	static double c0;									//�����g���i���������j


//*** [�R���X�g���N�^] ****************************
	public FeCr_PD_2D_001(){
		xwidth=400; yheight=400; 			//�`���ʂ̉��Əc�̒����i�s�N�Z���P�ʁj
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

//*** [main �v���O����] *****************************
	public static void main(String[] args) throws Exception{//��O�����͍s��Ȃ�

		FeCr_PD_2D_001 prog=new FeCr_PD_2D_001();	//FeCr_PD_2D_001�̃C���X�^���Xprog�𐶐�

		int i, j; 										//����
		int ip, im, jp, jm; 								//�����ii+1, i-1, j+1, j-1�j
		double delt;									//���ԍ��݁i�������j
		double al;										//�v�Z�̈�̂P�ӂ̒���
		double time1max;							//�v�Z���ԁi�ő�J�E���g���j
		double [][] ck = new double[ND][ND];		//�g�U�|�e���V����
		double [][] ch2 = new double[ND][ND]; 	//�g�D���̔Z�x�f�|�^�\���z��
		double mu_chem, mu_surf, mu_str, mu_mag;			//�e�|�e���V����
		double c1, c2;								//�Z�x(Fe:1, Cr:2)
		double a01, atomNo;			//Fe(bcc)�̊i�q�萔�ƁA�P�ʖE���̌��q��
		double vm0;													//molar volume

		double L0;									//���q�ԑ��ݍ�p�p�����[�^
		double L0_0, L0_1;					//���q�ԑ��ݍ�p�p�����[�^���̌W��
		double kappa_c;								//�Z�x���z�G�l���M�[�W��

		double eta;										//�i�q�~�X�}�b�`
		double c11a, c12a, c44a;			//Fe(bcc)�̒e���萔
		double c11b, c12b, c44b;			//Cr(bcc)�̒e���萔
		double c11, c12, c44;					//�S�̂̒e���萔
		double y100;									//�e���G�l���M�[�֐�

		//���C�ߏ�G�l���M�|�֘A�p�����|�^
		double Tc, d2Tc;												//�L�����|���x�A����т��̑g���ɑ΂������
		double tau;														//�L�����|���x�ŋK�i���������x
		double ftau, dftau;											//���C�ߏ�G�l���M�|�֐��Ɖ��x�ɑ΂������
		double Bc, d2Bc;											//�{�|�A���q�A����т��̑g���ɑ΂������
		double TcFe, TcCr, TcCrFe0, TcCrFe1; //�L�����[���x�֘A�̌W��
		double BcFe, BcCr, BcCrFe0; 				//�{�[�A���q�֘A�̌W��
		double p_mag, D_mag;										//p��D

		double Mc;									//�Փ��x�֐��Ƃ��̔���
		double c_flu;									//�Z�x��̗h�炬�̑傫��
		double cddtt;									//�Z�x�̑���

		double b1;									//�����u���b�N�T�C�Y
		double c2ip, c2im, c2jp, c2jm; 			//�����u���b�N�ɂ�����c2�𒆐S�ɁA���̏㉺���E�̔Z�x
		double sumc, dc; 							//�Z�x��̑��a�A���ϑg������̂���


//---- �e��p�����[�^�ݒ� ----------------------------------------------------
		temp=673.0;  			//[K]
		c0=0.5;					//�����̕��ϑg���i�����ł�A-30at%B������ݒ肵�Ă���j
		delt=0.02;					//���Ԃ�����
		time1=0.0;				//�v�Z�̌J�Ԃ���
		time1max=1.0e+08;		//�J�Ԃ��񐔂̍ő�l

		al=30.0; 					//�Q�����v�Z�̈�̂P�ӂ̒���(nm)
		al=al*1.0e-9;				//(m)�ɕϊ�
		b1=al/(double)nd;		//�����P�u���b�N�̃T�C�Y

		a01=0.28664e-09;			//�i�q�萔(bcc Fe) (m)
		atomNo=2.0;			//�P�ʖE���̌��q��
		vm0=6.02E23*a01*a01*a01/atomNo;//�����̐�

		L0_0=20500.0;		//���q�ԑ��ݍ�p�p�����[�^�iJ/mol�jFe-Cr(bcc)
		L0_1=-9.68;		//L=L0_0+L0_1*T
		L0=(L0_0+L0_1*temp)/RR/temp; 			//���q�ԑ��ݍ�p�p�����|�^�̖�������

		kappa_c=2.0e-15;			//�Z�x���z�G�l���M�|�W���A�P�ʂ�[J m^2/mol]
		kappa_c=kappa_c/b1/b1/RR/temp;		//�Z�x���z�G�l���M�|�W���A(b1^2*rr*temp)�Ŗ�������

		eta=0.00614;	//�i�q�~�X�}�b�`

		c11a=2.3310e+11; 	//bccFe�̒e���萔(Pa)
		c12a=1.3544e+11;
		c44a=1.1783e+11;

		c11b=3.500e+11; 		//bccCr�̒e���萔(Pa)
		c12b=0.678e+11;
		c44b=1.008e+11;

		c11a=c11a*vm0/RR/temp; 						//RT�Ŗ�������
		c12a=c12a*vm0/RR/temp;
		c44a=c44a*vm0/RR/temp;
		c11b=c11b*vm0/RR/temp;
		c12b=c12b*vm0/RR/temp;
		c44b=c44b*vm0/RR/temp;

		c11=(1.0-c0)*c11a+c0*c11b;					//�����̒e����
		c12=(1.0-c0)*c12a+c0*c12b;
		c44=(1.0-c0)*c44a+c0*c44b;

		y100=c11+c12-2.0*(c12*c12/c11);			//�e�����̊֐�Y<100>

		TcFe=1043.0; TcCr=-311.5; TcCrFe0=1650.0; TcCrFe1=550.0;//�L�����[���x�֘A�̌W��
		BcFe=2.22;   BcCr=-0.008; BcCrFe0=-0.85;//�{�[�A���q�֘A�̌W��
		p_mag=0.4;											//(bcc)
		D_mag=518.0/1125.0+11692.0/15975.0*(1.0/p_mag-1.0);

		Mc=c0*(1.0-c0);			//�g�U�̈Փ��x
		c_flu=0.1;

//---- ����0�ɂ����鏉���Z�x��ݒ� ----------------------------------------------------
		prog.ini_comp_field();
		//prog.datin();				//�t�@�C������ǂݍ��ޏꍇ

//---- �Z�x��̎��Ԕ��W�̌v�Z ----------------------------------------------------
		while(time1<=time1max){

//---- �Z�x��̕\�� ----------------------------------------------------
			//�J�E���g����200�̔{�������ɔZ�x��`��
			if((((int)(time1) % 100)==0)){ prog.update_draw(g); }	//�`�掞�̃`���c�L��}���邽��
			//if((((int)(time1) % 200)==0)){ prog.repaint(); }

//---- �Z�x��̕ۑ� ----------------------------------------------------
			if((((int)(time1) % 200)==0)){ prog.datsave(); }			//�J�E���g����200�̔{�������ɔZ�x��ۑ�
			//if(time1==3000.0){ prog.datsave(); }						//�J�E���g����3000�̎��ɔZ�x��ۑ�

//---- �g�U�|�e���V�����̌v�Z -----------------------------------------------------
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1; jp=j+1; jm=j-1;
					if(i==ndm){ip=0;}  if(i==0){im=ndm;} 	//�����I���E����
					if(j==ndm){jp=0;}  if(j==0){jm=ndm;} 	//�����I���E����

					c2=ch[i][j]; 		c1=1.0-c2; 				//�ʒu(i,j)�ɂ�����c1��c2
					c2ip=ch[ip][j]; c2im=ch[im][j]; 
					c2jp=ch[i][jp]; c2jm=ch[i][jm]; //�����ɂ�����c2�̑O�㍶�E�̔Z�x

		 			mu_chem=L0*(c1-c2)+Math.log(c2)-Math.log(c1); 				//���w�|�e���V������
					mu_surf=-2.0*kappa_c*(c2ip+c2im+c2jp+c2jm-4.0*c2);		//�Z�x���z�̃|�e���V����
					mu_str=2.0*eta*eta*y100*(c2-c0); 								//�e���|�e���V����

					//Tc=TcFe*c1+TcCr*c2;//�L�����[���x
					//d2Tc=-TcFe+TcCr;//Tc��Cr�g������
					//Bc=BcFe*c1+BcCr*c2;//�P���q������̎����̋����i�{�[�A���q�Ŗ��������j
					//d2Bc=-BcFe+BcCr;//Bc��Cr�g������

					Tc=TcFe*c1+TcCr*c2+c1*c2*(TcCrFe0+TcCrFe1*(c2-c1));//�L�����[���x
					d2Tc=-TcFe+TcCr+2.0*TcCrFe1*c1*c2+(c1-c2)*(TcCrFe0+TcCrFe1*(c2-c1));//Tc��Cr�g������
					Bc=BcFe*c1+BcCr*c2+BcCrFe0*c1*c2;//�P���q������̎����̋����i�{�[�A���q�Ŗ��������j
					d2Bc=-BcFe+BcCr+BcCrFe0*(c1-c2);//Bc��Cr�g������
					if(Tc<0.0){Tc=-Tc;  d2Tc=-d2Tc;}
					if(Bc<0.0){Bc=-Bc;  d2Bc=-d2Bc;}

					tau=temp/Tc; 															//�т̒�`
					if(tau<=1.0){ 
						//f(��)�ƁAf(��)�̃є����̌v�Z
						ftau=1.0-1.0/D_mag*(79.0/140.0/p_mag/tau
															+474.0/497.0*(1.0/p_mag-1.0)*( Math.pow(tau,3.0)/6.0
															+Math.pow(tau,9.0)/135.0+Math.pow(tau,15.0)/600.0) );
						dftau=-1.0/D_mag*(-79.0/140.0/p_mag/tau/tau
															+474.0/497.0*(1.0/p_mag-1.0)*( Math.pow(tau,2.0)/2.0
															+Math.pow(tau,8.0)/15.0+Math.pow(tau,14.0)/40.0) );
					}
					else{	
						ftau=-1.0/D_mag*(Math.pow(tau,-5.0)/10.0+Math.pow(tau,-15.0)/315.0
														+Math.pow(tau,-25.0)/1500.0);
						dftau=1.0/D_mag*(Math.pow(tau,-6.0)/2.0+Math.pow(tau,-16.0)/21.0
														+Math.pow(tau,-26.0)/60.0);
					}

					mu_mag=ftau/(Bc+1.0)*d2Bc-dftau*tau/Tc*d2Tc*Math.log(Bc+1.0); 
																							//���C�ߏ�G�l���M�[�̃|�e���V����

					ck[i][j]=mu_chem+mu_surf+mu_str+mu_mag;			//�g�U�|�e���V����

				}
			}

//---- �Z�x��̎��ԕω�(����`�g�U�������̍�����z��@) ------------------------------------
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
					if(i==ndm) {ip=0;}	if(i==0) {im=ndm;} 	//�����I���E����
					if(j==ndm) {jp=0;}	if(j==0) {jm=ndm;} 	//�����I���E����
					cddtt=Mc*(ck[ip][j]+ck[im][j]+ck[i][jp]+ck[i][jm]-4.0*ck[i][j]);		//����`�g�U������
					//ch2[i][j]=ch[i][j]+cddtt*delt; 										//�Z�x��̎��Ԕ��W
					ch2[i][j]=ch[i][j]+( cddtt+c_flu*(2.0*Math.random()-1.0) ) *delt;	//�Z�x��̎��Ԕ��W�i�Z�x�h�炬�����j
				}
			}

//*** [�Z�x��̎��x�̕␳] *******************************************************
//*** ���l�v�Z�ł���̂ŔZ�x��̎��x�̕␳���s���i���ۂɂ͖��X�e�b�v�s���K�v�͂Ȃ��j�B] ****
  		sumc=0.;
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					sumc=sumc+ch2[i][j]; 					//�Z�x��̐ϕ�
				}
			}
			dc=sumc/(double)nd/(double)nd-c0;			//�Z�x��̕ϓ���

			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ch[i][j]=ch2[i][j]-dc; 						//�Z�x��̕␳
					if(ch[i][j]>=1.){ch[i][j]=1.0-1.0e-6;} 	//�Z�x��1�𒴂����ꍇ�̕␳
					if(ch[i][j]<=0.){ch[i][j]=1.0e-6;} 			//�Z�x��0��؂����ꍇ�̕␳
				}
			}

//******[���ԑ���]*************************************************
			time1=time1+1.0; 	//�v�Z���Ԃ̑���
		}//while

//----------------------------------------------------------------

		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");

}//main

// �ȉ��̓T�u���[�`���ł���B
// **** �����Z�x��ݒ� *****************************************************
	public void ini_comp_field(){
		int i, j;
		double rnd0, fac1;

		fac1=0.01; 	//�����Z�x�h�炬�̍ő�ω��ʂ�1%�ɐݒ�
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				rnd0=2.0*Math.random()-1.0;  ch[i][j]=c0+rnd0*fac1; 	//�����ɂď����Z�x���ݒ�
			}
		}

	}

// *********************************************************
//�`���c�L�h�~�̂���draw�֐���ʂɒ�`
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** �Z�x��`�� *****************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);		//Window���N���A

		int i, j, ii, jj;
		double c, x, xmax, xmin, y, ymax, ymin, rad0;
		int ixmin=0, iymin=0, igx, igy, irad0;
		int ixmax=xwidth, iymax=yheight;
		int icol;

		xmin=0.; xmax=1.; 						//�����̍ŏ��l�A�ő�l
		ymin=0.; ymax=1.; 						//�c���̍ŏ��l�A�ő�l

		rad0=1.0/(double)nd/2.0; 				//�����u���b�N�̒����̔���
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 ); 	//rad0�̃s�N�Z����

		System.out.printf("%f \n", time1); 	//�v�Z�̌J�Ԃ��񐔂�W�����o�͂ɕ\��

		for(i=0;i<=nd;i++){
			for(j=0;j<=nd;j++){
				//�Z�x��̈ʒu���W�i���ۂ̒l�j
				x=1.0/(double)nd*(double)i+rad0;
				y=1.0/(double)nd*(double)j+rad0;
				//�Z�x��̈ʒu���W�i�X�N���[�����W�ɕϊ��j
				igx=(int)( ((double)ixmax-(double)ixmin)*(x-xmin)/(xmax-xmin)+(double)ixmin );
				igy=(int)( ((double)iymax-(double)iymin)*(y-ymin)/(ymax-ymin)+(double)iymin );

				//�X�̍����u���b�N�̔Z�x�l
				ii=i; jj=j;
				if(i==nd){ii=0;} if(j==nd){jj=0;}											//�����I���E����
				icol=(int)(255.0*(1.0-ch[ii][jj]));											//�F���~�����O���[�X�P�[���ɂ���
				//icol=(int)(255.0*ch[ii][jj]);												//���Â𔽓]����ꍇ
				if(icol>=255){icol=255;} if(icol<=0){icol=0;}							//���Â͈̔͂̕␳
				g.setColor(new Color(icol,icol,icol)); 									//�Z�x�𖾈ÂŐݒ�
				g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);		//�X�̍����u���b�N�`��
			}
		}
	}

//*** [�f�|�^�̕ۑ�] ************************************
	private void datsave() throws Exception{
		int	i, j;

		//�ۑ��t�@�C������test.dat�Ƃ���B
		PrintWriter outfile= new PrintWriter(
			new BufferedWriter(new FileWriter("test.dat", true)) );		//�t�@�C���̃I�[�v���ǋL

		outfile.println(time1);				//�J�E���g�̏�������
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				outfile.println(ch[i][j]);		//�Z�x��̏�������
			}
		}
		outfile.close();						//�t�@�C���̃N���[�Y
	}

//*** [�f�|�^�̓Ǎ���] ************************************
	private void datin() throws Exception{
		int	i, j;
		String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//�t�@�C���̃I�[�v��

		s_data=infile.readLine();  									//����������ēǂݍ���
		time1=new Double(s_data).doubleValue();				//�����𐔒l�֕ϊ�
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				s_data=infile.readLine();  							//����������ēǂݍ���
				ch[i][j]=new Double(s_data).doubleValue();		//�����𐔒l�֕ϊ�
			}
		}
		infile.close();													//�t�@�C���̃N���[�Y
	}

//****************************************************************
}//FeCr_PD_2D_001
//*** �v���O�����I�� ************************************************************
