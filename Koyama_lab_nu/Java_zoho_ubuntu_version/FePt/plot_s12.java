import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class plot_s12 extends Frame{

	static int ND=128;					//�g�D�P�ӂ̕�����
	static int nd=ND, ndm=ND-1;
	static int width;					// Window�S�̂̕�
	static int height;					// Window�S�̂̍���
	static int xwidth;					// �`��̈�̕�
	static int yheight;					// �`��̈�̍���
	static int insetx;					// Window�̘g�̕��i���E����щ��j
	static int insety;					// Window�̘g�̕��i��j
	static double [][] s1h=new double[ND][ND];		//�g�D���̒����ϐ��f�|�^s1�z��
	static double [][] s2h=new double[ND][ND];		//�g�D���̒����ϐ��f�|�^s2�z��

	//�O���[�o���ϐ���
	static double time1;	//�v�Z���ԁi�J�E���g���j
	Image buff;
	static Graphics bg, g;

	//*************** �R���X�g���N�^ ****************************
	public plot_s12(){
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

//************** main *****************************
	public static void main(String[] args) throws Exception{//��O�����͍s��Ȃ�

		plot_s12 prog=new plot_s12();

		int  i, j, k, l; 						//����
		double sumc;
    String str_data, str_time1;

		BufferedReader infile=new BufferedReader(new FileReader("test_FePt.dat"));
		while( (str_time1=infile.readLine())!=null ){
			time1=new Double(str_time1).doubleValue();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					str_data=infile.readLine();
					s1h[i][j]=new Double(str_data).doubleValue();
					str_data=infile.readLine();
					s2h[i][j]=new Double(str_data).doubleValue();
				}
			}
			prog.update_draw(g);
			//if(time1==10000.0){ prog.datsave(); }		//�J�E���g�����w�肵�ď��ۑ�
		}

		infile.close();

//----------------------------------------------------------------
		System.out.printf("\n �I�����܂����B�O���t�̉E��~���N���b�N���ďI�����Ă��������B\n");

	}//main

// *******************************************************************************
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

				icol_r=(int)(255.*c_r);  	icol_g=(int)(255.*c_g);  icol_b=(int)(255.*c_b);		//256�K�w�ɕϊ�
				g.setColor(new Color(icol_r, icol_g, icol_b)); 									//�����u���b�N�̐F��ݒ�
				g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);					//�X�̍����u���b�N�`��
			}
		}
	}

//**** [�f�|�^�̕ۑ�] ************************************
	private void datsave() throws Exception{
		int	i, j;

		//�ۑ��t�@�C������test.dat�Ƃ���B
		PrintWriter outfile= new PrintWriter(							//�t�@�C���̃I�[�v��
			new BufferedWriter(new FileWriter("test_selected.dat", true)) );	//�ǋL

		outfile.println(time1);											//�J�E���g�̏�������
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//outfile.printf("%e  %e  ",s1h[i][j], s2h[i][j]);									//phase field�̏�������
				outfile.println(s1h[i][j]);									//phase field�̏�������
				outfile.println(s2h[i][j]);									//phase field�̏�������
			}
		}
		outfile.close();													//�t�@�C���̃N���[�Y
	}

//****************************************************************
}