//*** [�v���O���� (MT2D_001.java)] ************************************************
//*** [�C���|�[�g��] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class MT2D_001 extends Frame{
//*** [�O���[�o���ϐ�] *****************************************************************************************
	static int ND=128;				//�g�D�P�ӂ̕�����
	static int IG=7;					//2^IG=ND
	static int nd=ND;					//�v�Z�̈�̂P�ӂ̕�����
	static int ndm=ND-1;			//�v�Z�̈�̂P�ӂ̕�����-1
	static int nd2=ND/2;				//�v�Z�̈�̂P�ӂ̕�����/2
	static int ig=IG;					//2^ig=ND
	static int width;					// Window�S�̂̕�
	static int height;					// Window�S�̂̍���
	static int xwidth;					// �`��̈�̕�
	static int yheight;					// �`��̈�̍���
	static int insetx;					// Window�̘g�̕��i���E����щ��j
	static int insety;					// Window�̘g�̕��i��j
	static double PI=3.141592;		//��
	static double RR=8.3145;		//�K�X�萔
	static double [][] s1h=new double[ND][ND];		//�g�D���̒����ϐ��f�|�^s1�z��
	static double [][] s2h=new double[ND][ND];		//�g�D���̒����ϐ��f�|�^s2�z��
	static Graphics g;									//���R�G�l���M�[�Ȑ���ʂ̃O���t�B�b�N�X�I�u�W�F�N�g
	static double time1;									//�v�Z���ԁi�J�E���g���j
	static double temp; 									//���x[K]

	//�����t�[���G�ϊ��̃v���O�����ɂ��ẮA����(9)�����Q�Ɖ������B
	static double qs;										//�t�|���G�ϊ�(qs:-1)�Ƌt�t�|���G�ϊ�(qs:1)�̋��
	static double [][] xr = new double[ND][ND];		//�t�|���G�ϊ��̎����p�|�g�Ɏg�p����z��
	static double [][] xi = new double[ND][ND];		//�t�|���G�ϊ��̋����p�|�g�Ɏg�p����z��
	static double [] xrf = new double[ND];				//�t�|���G�ϊ��̎����p�|�g�Ɏg�p����z��
	static double [] xif = new double[ND];				//�t�|���G�ϊ��̋����p�|�g�Ɏg�p����z��
	static double [] s = new double[ND];				//sin�̃e�[�u��
	static double [] c = new double[ND];				//cos�̃e�[�u��
	static int [] ik = new int[ND];						//�r�b�g���]����̔z��

//*** [�R���X�g���N�^] ****************************
	public MT2D_001(){
		xwidth=400; yheight=400; 			//�`���ʂ̉��Əc�̒����i�s�N�Z���P�ʁj
		insetx=4; insety=30;					//�`���ʂ̂ӂ��̒���
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

		MT2D_001 prog=new MT2D_001();		//MT2D_001�̃C���X�^���Xprog�𐶐�

		double [][] ec11 = new double[ND][ND];			//�S���c�ϓ���
		double [][] ec22 = new double[ND][ND];			//�S���c�ϓ���
		double [][] ep11h0 = new double[ND][ND];		//�ϑԘc
		double [][] ep22h0 = new double[ND][ND];		//�ϑԘc
		double [][] ep11qrh0 = new double[ND][ND];		//�S���c�ϓ��ʂ̃t�[���G�ϊ��i�������j
		double [][] ep11qih0 = new double[ND][ND];		//�S���c�ϓ��ʂ̃t�[���G�ϊ��i�������j
		double [][] ep22qrh0 = new double[ND][ND];		//�S���c�ϓ��ʂ̃t�[���G�ϊ��i�������j
		double [][] ep22qih0 = new double[ND][ND];		//�S���c�ϓ��ʂ̃t�[���G�ϊ��i�������j
		double [][] eta_s1 = new double[4][4];				//�o���A���g�P�̕ϑԘc
		double [][] eta_s2 = new double[4][4];				//�o���A���g�Q�̕ϑԘc
		double [][] s1k_su = new double[ND][ND];		//���z�|�e���V����
		double [][] s2k_su = new double[ND][ND];		//���z�|�e���V����

		double s1, s2;										//�}���e���T�C�g��phase field
		double s1k_chem, s1k_str;							//���w�|�e���V�����A�e���|�e���V����
		double s2k_chem, s2k_str;							//���w�|�e���V�����A�e���|�e���V����
		double c11, c12, c44, lam0, mu0, nu0;			//�e���萔
		double el_fac;										//�e���萔�K�i���ϐ�
		double ep11T, ep22T;								//�X�̘c�����̘a
		double ep11_0, ep22_0;							//�g�D���̕ϑԘc�̕��ϒl
		double ep11_a, ep22_a, ep12_a, ep21_a;		//�O�͂ɋN������c
		double sig11_a, sig22_a;							//�O��
		double Z11ep, Z12ep, Z21ep, Z22ep;			//�t�[���G�ϊ����Ɏg�p����W��
		double sum11, sum22;
		double s1ddtt, s2ddtt;								//phase field�̎��ԕω���
		double delt;											//���Ԃ�����

		int   i, j, k, l, ii=0, jj=0, kk, iii, jjj;								//����
		int   p, q, m, n;										//����
		int   ip, im, jp, jm;									//����
		double al, temp;										//�v�Z�̈�A���x
		double time1max;									//�ő厞�ԁi�v�Z���~�߂�ۂɎg�p�j
		double b1, vm0, atom_n;							//�K�i�������A�����̐ρA�P�ʖE���̌��q��
		double smob;										//�}���e���T�C�g�ϑԃ_�C�i�~�N�X�̊ɘa�W��
		double nx, ny, nxx, nyy, alnn;						//�t��Ԃ̊�{�x�N�g���A���̎���A�m����

		double AA0, AA1, AA2, AA3;						//���w�I�쓮�͒萔
		double a1_c, b1_c, c1_c;							//�i�q�萔
		double kappa_s1, kappa_s2;						//���z�G�l���M�|�萔
		double ds_fac;										//phase field�̗h�炬�̑傫��

//---- �e��p�����[�^�ݒ� ----------------------------------------------------

		delt=0.1;	//���Ԃ����ݓ���

		temp=500.0;						//���x[K]
		al=250.0*1.0E-09;				//�v�Z�̈�[m]
		b1=al/nd;							//�����u���b�N�̒���[m]

		time1=-10.0;						//�����ݒ莞��
		time1max=1.0+1.0e+07;		//�v�Z���Ԃ̍ő�l

		smob=1.0;						//�}���e���T�C�g�ϑԃ_�C�i�~�N�X�̊ɘa�W��
		ds_fac=0.01;						//phase field�̗h�炬�W��

		AA0=1000.0;						//�}���e���T�C�g�ϑԂ̉��w�I�쓮��[J/mol]
		AA0=AA0/RR/temp;				//��������
		AA1=10.0; AA2=3.0*AA1+12.0; AA3=2.0*AA1+12.0;	//���w�I�쓮�͒萔

		kappa_s1=5.0e-15;												//���z�G�l���M�|�萔[Jm^2/mol]
		kappa_s1=kappa_s2=kappa_s1/RR/temp/b1/b1;			//��������

		a1_c=b1_c=c1_c=3.563e-10;									//�i�q�萔[nm]
		atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//�����̐ς̌v�Z�ifcc������j[m^3/mol]

//--- s1��̕ϑԘc�̐ݒ� ---
		eta_s1[1][1]=0.08; eta_s1[2][2]=-0.04; eta_s1[3][3]=0.0;
		eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.0;

//--- s2��̕ϑԘc�̐ݒ� ---
		eta_s2[1][1]=eta_s1[2][2]; 	eta_s2[2][2]=eta_s1[1][1]; 	eta_s2[3][3]=0.0;
		eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.0;

//---�e���萔(Ni�̏ꍇ) ------------------------------------
		el_fac=1.0E+11*vm0/RR/temp;
		c11=2.508*el_fac;
		c44=1.235*el_fac;
		c12=1.500*el_fac;
		//c12=c11-2.0*c44;
		lam0=c12;	mu0=c44;			//���[���̒萔
		nu0=lam0/2.0/(lam0+mu0);	//�|�A�\����

//---- �O�͂̐ݒ� ---------------------------------------------
		sig22_a=-200.0*1.0e+06*vm0/RR/temp;						//�����ł͊O�͂�0�ɐݒ�
		ep11_a=-lam0/4.0/mu0/(lam0+mu0)*sig22_a;					//���ʘc��z��
		ep22_a=(lam0+2.0*mu0)/4.0/mu0/(lam0+mu0)*sig22_a;
		ep12_a=ep21_a=0.0;

//**** phase field�̏�����ݒ�ƁAsin�����cos�e�[�u���̐ݒ� *******************************************
		prog.ini_field();				//phase field�̏�����ݒ�
		prog.table();					//�t�[���G�ϊ��̂��߂�sin��cos�̃e�[�u���ƃr�b�g�ϊ��z��̐ݒ�

//**** phase field�̎��Ԕ��W�̌v�Z *******************************************************************************
		while(time1<=time1max){

//---- phase field�̕\�� ----------------------------------------------------
			if((((int)(time1) % 50)==0)){ prog.update_draw(g); }		//�J�E���g����50�̔{�������ɏ��`��
			//if((((int)(time1) % 100)==0)){ prog.repaint(); } 

//---- phase field�̕ۑ� ----------------------------------------------------
			//if((((int)(time1) % 200)==0)){ prog.datsave(); }			//�J�E���g����200�̔{�������ɏ��ۑ�
			//if(time1==3000.0){ prog.datsave(); }						//�J�E���g����3000�̎��ɏ��ۑ�

//***** ���z�|�e���V���� *******************************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1; jp=j+1; jm=j-1;
					if(i==ndm){ip=0;}  if(i==0){im=ndm;}
					if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
					s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1h[i][j]);		//��(4)
					s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]-4.0*s2h[i][j]);		//��(4)
				}
			}

//**** �ϑԘc��̃t�|���G�ϊ� ep11 ********************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]+eta_s2[1][1]*s2h[i][j];		//��(5)
					xi[i][j]=0.0;
				}
			}
			qs=-1.0; prog.rcfft();		//����Ԃ���t�[���G��Ԃւ̕ϊ�(qs<0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ 
					ep11qrh0[i][j]=xr[i][j];  ep11qih0[i][j]=xi[i][j];
				}
			}
			ep11qrh0[0][0]=ep11qih0[0][0]=0.0;

//**** �ϑԘc��̃t�|���G�ϊ� ep22 ********************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]+eta_s2[2][2]*s2h[i][j]; 	//��(5)
					xi[i][j]=0.0;
				}
			}
			qs=-1.0; prog.rcfft();		//����Ԃ���t�[���G��Ԃւ̕ϊ�(qs<0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ 
					ep22qrh0[i][j]=xr[i][j];  ep22qih0[i][j]=xi[i][j];
				}
			}
			ep22qrh0[0][0]=ep22qih0[0][0]=0.0;

//*** �ϑԘc��̕��ϒl�̎Z�o ***
			sum11=sum22=0.0;
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j]; }
			}
			ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** �S���c�ϓ���ec11�̌v�Z *************************************
			for(i=0;i<=ndm;i++){
				if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
				for(j=0;j<=ndm;j++){
					if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
					alnn=Math.sqrt((double)ii*(double)ii+(double)jj*(double)jj);
					if(alnn==0.){alnn=1.;}
					nxx=(double)ii/alnn*(double)ii/alnn;
					nyy=(double)jj/alnn*(double)jj/alnn;
					Z11ep=nxx*(2.0*(1.0-nu0)-nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
					Z12ep=nxx*(2.0*nu0     -nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
					xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j]; 		//��(10)
					xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j]; 		//��(10)
				}
			}
			qs=1.0; prog.rcfft();		//�t�[���G��Ԃ������Ԃւ̕ϊ�(qs>0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ ec11[i][j]=xr[i][j]; }
			}

//***** �S���c�ϓ���ec22�̌v�Z *****************************
			for(i=0;i<=ndm;i++){
				if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
				for(j=0;j<=ndm;j++){
					if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
					alnn=Math.sqrt((double)ii*(double)ii+(double)jj*(double)jj);
					if(alnn==0.){alnn=1.;}
					nxx=(double)ii/alnn*(double)ii/alnn;
					nyy=(double)jj/alnn*(double)jj/alnn;
					Z21ep=nyy*(2.0*nu0     -nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
					Z22ep=nyy*(2.0*(1.0-nu0)-nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
					xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j]; 		//��(10)
					xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j]; 		//��(10)
				}
			}
			qs=1.0; prog.rcfft();		//�t�[���G��Ԃ������Ԃւ̕ϊ�(qs>0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ ec22[i][j]=xr[i][j]; }
			}

//****** �|�e���V�����Ɣ��W�������̌v�Z ************************************************************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					s1=s1h[i][j];  	s2=s2h[i][j];

//****** ���w�|�e���V�����̌v�Z ********************************
					s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2)); 		//��(2)
					s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2)); 		//��(2)

//****** �e���|�e���V�����̌v�Z ********************************
					ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
					ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;

					s1k_str=ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
								 +ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]); 		//��(7)
					s2k_str=ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
								 +ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]); 		//��(7)

//****** phase field�̎��Ԕ��W�̌v�Z ********************************
					s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str); 				//��(12)
					s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str); 				//��(12)
					s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*Math.random()-1.0) )*delt;
					s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*Math.random()-1.0) )*delt;

//--- s�̕ψ�(0<=s<=1)�̕␳ ---
					if(s1h[i][j]>=1.0){s1h[i][j]=1.0;}  if(s1h[i][j]<=0.0){s1h[i][j]=0.0;}
					if(s2h[i][j]>=1.0){s2h[i][j]=1.0;}  if(s2h[i][j]<=0.0){s2h[i][j]=0.0;}
				}
			}

//**** [���ԑ���] *************************************************
			time1=time1+1.0;
		}//while

//----------------------------------------------------------------
		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");
	}//main

// �ȉ��̓T�u���[�`���ł���B
// **** [phase field�̏����ݒ�] *****************************************************
	public void ini_field(){
		int i, j;
		double fac1;

		fac1=0.5; //�ő�̏����h�炬
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//s1h[i][j]=fac1*Math.random(); s2h[i][j]=fac1*Math.random();						//�ψ�Ɋj��u���ꍇ
				if(Math.abs(j-nd2)<(nd/40)){s1h[i][j]=Math.random(); s2h[i][j]=Math.random();}	//�����Ɋj��u���ꍇ
			}
		}

	}

	//******* [�t�[���G�ϊ��̂��߂�sin��cos�̃e�[�u���ƁA�r�b�g�ϊ��z��̐ݒ�] **********************
	public void table(){
		int it, it1, it2, mc, mn;
		double q;

		q=2.0*PI/(double)nd;
		for(it=0;it<=nd2-1;it++){ c[it]=Math.cos(q*(double)it); s[it]=Math.sin(q*(double)it); }

		ik[0]=0; mn=nd2; mc=1;
		for(it1=1;it1<=ig;it1++){
			for(it2=0;it2<=mc-1;it2++){ ik[it2+mc]=ik[it2]+mn; }
			mn=mn/2; mc=2*mc;
		}
	}

// *******************************************************************************
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** [phase field�̕`��] *******************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);//Window���N���A

		int i, j, ii, jj;
		int icol, icol_r, icol_g, icol_b;
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

				c_r= s1h[ii][jj];							//s1���
				c_g= s2h[ii][jj]; 							//s2���
				c_b=1.0-c_r-c_g; 						//�ϑԑO�̑����
				if(c_r>1.0){c_r=1.0;}  if(c_r<0.0){c_r=0.0;}
				if(c_g>1.0){c_g=1.0;}  if(c_g<0.0){c_g=0.0;}
				if(c_b>1.0){c_b=1.0;}  if(c_b<0.0){c_b=0.0;}

				icol_r=(int)(255.*c_r);  	icol_g=(int)(255.*c_g);  icol_b=(int)(255.*c_b);		//256�K�w�ɕϊ�
				g.setColor(new Color(icol_r, icol_g, icol_b)); 									//�����u���b�N�̐F��ݒ�
				g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);					//�X�̍����u���b�N�`��
			}
		}
	}

//***** [�P���������t�[���G�ϊ�(FFT)] **************************************
	public void fft(){
		int ix, ka, kb, l2, lf, mf, n2, nf;
		double tj, tr;

		l2=1;
		for(lf=1;lf<=ig;lf++){
			n2=nd2/l2;
			for(mf=1;mf<=l2;mf++){
				for(nf=0;nf<=n2-1;nf++){
					ix=nf*l2; ka=nf+2*n2*(mf-1); kb=ka+n2;
					tr=xrf[ka]-xrf[kb];            tj=xif[ka]-xif[kb];
					xrf[ka]=xrf[ka]+xrf[kb];       xif[ka]=xif[ka]+xif[kb];
					xrf[kb]=tr*c[ix]+tj*qs*s[ix];    xif[kb]=tj*c[ix]-tr*qs*s[ix];
				}
			}
			l2=l2*2;
		}
	}

//**** [�Q���������t�[���G�ϊ�(RC FFT)] ***********************************
	public void rcfft(){
		int i, ic, ir, j;

		for(ir=0;ir<=ndm;ir++){
			for(ic=0;ic<=ndm;ic++){ xrf[ic]=xr[ir][ic];  xif[ic]=xi[ir][ic]; }
			fft();
			for(ic=0;ic<=ndm;ic++){ xr[ir][ic]=xrf[ik[ic]];  xi[ir][ic]=xif[ik[ic]]; }
		}
		for(ic=0;ic<=ndm;ic++){
			for(ir=0;ir<=ndm;ir++){ xrf[ir]=xr[ir][ic];  xif[ir]=xi[ir][ic]; }
			fft();
			for(ir=0;ir<=ndm;ir++){ xr[ir][ic]=xrf[ik[ir]];  xi[ir][ic]=xif[ik[ir]]; }
		}
		if(qs>0.0){return;}
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){ xr[i][j]=xr[i][j]/nd/nd;  xi[i][j]=xi[i][j]/nd/nd; }
		}
	}

//**** [�f�|�^�̕ۑ�] ************************************
	private void datsave() throws Exception{
		int	i, j;

		//�ۑ��t�@�C������test.dat�Ƃ���B
		PrintWriter outfile= new PrintWriter(							//�t�@�C���̃I�[�v��
			new BufferedWriter(new FileWriter("test.dat", true)) );	//�ǋL

		outfile.println(time1);											//�J�E���g�̏�������
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				outfile.println(s1h[i][j]);									//phase field�̏�������
				outfile.println(s2h[i][j]);									//phase field�̏�������
			}
		}
		outfile.close();													//�t�@�C���̃N���[�Y
	}

//**** [�f�|�^�̓Ǎ���] ************************************
	private void datin() throws Exception{
		int	i, j;
		String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//�t�@�C���̃I�[�v��

		s_data=infile.readLine();  										//������Ƃ��ēǂݍ���
		time1=new Double(s_data).doubleValue();					//�����𐔒l�֕ϊ�
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				s_data=infile.readLine();  								//������Ƃ��ēǂݍ���
				s1h[i][j]=new Double(s_data).doubleValue();			//�����𐔒l�֕ϊ�
				s_data=infile.readLine(); 								//������Ƃ��ēǂݍ���
				s2h[i][j]=new Double(s_data).doubleValue();			//�����𐔒l�֕ϊ�
			}
		}
		infile.close();														//�t�@�C���̃N���[�Y
	}

//*******************************************************************************************************
}//MT2D_001
//*** �v���O�����I�� ********************************************************************************
