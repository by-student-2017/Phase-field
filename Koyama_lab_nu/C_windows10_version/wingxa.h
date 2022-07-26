/* Functions for Gx using Win32API			mmvaviimxxx,mmviaximxii */

#ifndef	_WINGXA_H
#define	_WINGXA_H

#define		_CRT_SECURE_NO_WARNINGS	// printf等の警告を表示しない

#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<windows.h>
#include	<tchar.h>	// _T()
#include	<limits.h>

#define	m_pi	3.1415926535897932384626433832795					// π
#define	m_rad	0.017453292519943295769236907684886				// π/180

#define	_int(x)	((int)floor((x)+0.5))

static int	_wide=767,_hite=767,_xorg=383,_yorg=383,
			_leftmrg=1024-767,_topmrg=0;
static int	_red=255, _grn=255, _blu=255;
static int	_penwidth=0;
static double	_cx, _cy;	// 現在の位置(初期値は原点)

static HDC	_hdc,_hdc2;
static PAINTSTRUCT	_ghpaint;
static int	_sort,_close=0,_end=0;					// _end:描画終了フラグ
static HBITMAP	_hbmp;

static int	randi(int r){	return(rand() % r);}			// 整数の乱数

//////////////	Drawing
static void	gcolor(int r,int g,int b){	_red=r;	_grn=g; _blu=b;}

static void	glinewidth(int	w){	_penwidth=w;}

static void	gpnt(double x,double y){
	SetPixel(_hdc,_int(x)+_xorg,_hite-_int(y)-_yorg,RGB(_red,_grn,_blu));
}

static void	getpntcol(double x,double y,int *r,int *g,int *b){
	COLORREF	C=GetPixel(_hdc,_int(x)+_xorg,_hite-_int(y)-_yorg);
	*r=GetRValue(C); *g=GetGValue(C); *b=GetBValue(C);
}

static void	gline(double x,double y,double x1,double y1){
	HPEN	hpen=CreatePen(PS_SOLID,_penwidth,RGB(_red,_grn,_blu)),
			hpen_old=(HPEN)SelectObject(_hdc,hpen);
	POINT	p[2];
	p[0].x= _int(x )+_xorg;	p[0].y= _hite-_int(y )-_yorg;
	p[1].x= _int(x1)+_xorg;	p[1].y= _hite-_int(y1)-_yorg;
	Polyline(_hdc,p,2);
	SelectObject(_hdc,hpen_old);	DeleteObject(hpen);
}

static void	gmoveto(double x, double y){
	_cx = x;	_cy = y;
}

static void	glineto(double x, double y){
	gline(_cx, _cy, x, y);
	_cx = x;	_cy = y;
}

static void	draw_grid(double  d){
	if(d){
		int	i, max=(int)(_yorg/d);
		gcolor(0,0,100);
		for(i=-max; i<=max; i++){
			gline(-_xorg,i*d,_xorg,i*d);
			gline(i*d,-_yorg,i*d,_yorg);
		}
	}
	gcolor(100,100,100);
	gline(-_xorg,0,_xorg,0);
	gline(0,-_yorg,0,_yorg);
	gcolor(255,255,255);
}

static void	_gPolygon(POINT *p,int p_nmb,HBRUSH hbrush){
	HPEN	hpen=CreatePen(PS_SOLID,_penwidth,RGB(_red,_grn,_blu)),
			hpen_old=(HPEN)SelectObject(_hdc,hpen);
	HBRUSH	hbrush_old=(HBRUSH)SelectObject(_hdc,hbrush);
	Polygon(_hdc,p,p_nmb);
	SelectObject(_hdc,hbrush_old);	DeleteObject(hbrush);
	SelectObject(_hdc,hpen_old);	DeleteObject(hpen);
}

static void	gtriangle(double x,double y,double x1,double y1,double x2,double y2){
	HBRUSH	hbrush=(HBRUSH)GetStockObject(NULL_BRUSH);		// not painted
	POINT	p[3];
	p[0].x= _int(x )+_xorg;	p[0].y= _hite-_int(y )-_yorg;
	p[1].x= _int(x1)+_xorg;	p[1].y= _hite-_int(y1)-_yorg;
	p[2].x= _int(x2)+_xorg;	p[2].y= _hite-_int(y2)-_yorg;
	_gPolygon(p,3,hbrush);
}

static void	gtri(double x0,double y0,double x1,double y1,double x2,double y2){
	HBRUSH	hbrush=CreateSolidBrush(RGB(_red,_grn,_blu));		// painted
	POINT	p[3];
	p[0].x= _int(x0)+_xorg;	p[0].y= _hite-_int(y0)-_yorg;
	p[1].x= _int(x1)+_xorg;	p[1].y= _hite-_int(y1)-_yorg;
	p[2].x= _int(x2)+_xorg;	p[2].y= _hite-_int(y2)-_yorg;
	_gPolygon(p,3,hbrush);
}

static void	gtetragon(double x0,double y0,double x1,double y1,
					double x2,double y2,double x3,double y3){	// not painted
	HBRUSH	hbrush=(HBRUSH)GetStockObject(NULL_BRUSH);
	POINT	q[4];
	q[0].x= _int(x0)+_xorg;	q[0].y= _hite-_int(y0)-_yorg;
	q[1].x= _int(x1)+_xorg;	q[1].y= _hite-_int(y1)-_yorg;
	q[2].x= _int(x2)+_xorg;	q[2].y= _hite-_int(y2)-_yorg;
	q[3].x= _int(x3)+_xorg;	q[3].y= _hite-_int(y3)-_yorg;
	_gPolygon(q,4,hbrush);
}

static void	gtetra(double x0,double y0,double x1,double y1,
			double x2,double y2,double x3,double y3){		// painted
	HBRUSH	hbrush=CreateSolidBrush(RGB(_red,_grn,_blu));
	POINT	q[4];
	q[0].x= _int(x0)+_xorg;	q[0].y= _hite-_int(y0)-_yorg;
	q[1].x= _int(x1)+_xorg;	q[1].y= _hite-_int(y1)-_yorg;
	q[2].x= _int(x2)+_xorg;	q[2].y= _hite-_int(y2)-_yorg;
	q[3].x= _int(x3)+_xorg;	q[3].y= _hite-_int(y3)-_yorg;
	_gPolygon(q,4,hbrush);
}

static void	grectangle(double x0,double y0,double x1,double y1){// not painted
	gtetragon(x0,y0,x1,y0,x1,y1,x0,y1);
}

static void	grect(double x,double y,double u,double v){			// painted
	gtetra(x,y,u,y,u,v,x,v);
}

static void	_ellipse(double x,double y,double rx,double ry,HBRUSH hbrush){
	HPEN	hpen=CreatePen(PS_SOLID,_penwidth,RGB(_red,_grn,_blu)),
			hpen_old=(HPEN)SelectObject(_hdc,hpen);
	HBRUSH	hbrush_old=(HBRUSH)SelectObject(_hdc,hbrush);
	x+= _xorg;	y= _hite-y-_yorg;
	if(GetGraphicsMode(_hdc)==GM_COMPATIBLE)		// if Window NT
			Ellipse(_hdc,_int(x-rx),_int(y-ry),_int(x+rx),_int(y+ry));
	else	Ellipse(_hdc,_int(x-rx),_int(y-ry),_int(x+rx+1),_int(y+ry+1));
	SelectObject(_hdc,hbrush_old);DeleteObject(hbrush);
	SelectObject(_hdc,hpen_old);  DeleteObject(hpen);
}

static void	gellipse(double x,double y,double rx,double ry){
	HBRUSH	hbrush=(HBRUSH)GetStockObject(NULL_BRUSH);
	_ellipse(x,y,rx,ry,hbrush);
}

static void	gellip(double x,double y,double rx,double ry){
	HBRUSH	hbrush=CreateSolidBrush(RGB(_red,_grn,_blu));		// painted
	_ellipse(x,y,rx,ry,hbrush);
}

#define	gcircle(x,y,r)	gellipse(x,y,r,r)					// not painted
#define	gcirc(x,y,r)	gellip(x,y,r,r)							// painted

static void	color16(int c){
	int	col[16][3]={{0,0,0},{0,0,128},{204,77,51},{128,0,128},
					{0,128,0},{0,128,128},{128,128,0},{128,128,128},
					{178,178,178},{0,0,255},{255,0,0},{255,0,255},
					{0,255,0},{0,255,255},{255,255,0},{255,255,255}};
	gcolor(col[c][0],col[c][1],col[c][2]);
}

static void	gcls(void){
	SelectClipRgn(_hdc,NULL);	PatBlt(_hdc,0,0,_wide,_hite,BLACKNESS);
}

static void	gtextout(int x,int y,const char *s){
	SetBkMode(_hdc,TRANSPARENT);
	SetTextColor(_hdc,RGB(_red,_grn,_blu));
	TextOut(_hdc,x+_xorg,_hite-y-_yorg,(LPTSTR)s,(int)strlen(s));
}

//////////	other functions
static int	keypress(void){
	MSG	msg;
	if(PeekMessage(&msg,NULL,0,0,PM_REMOVE)){
		TranslateMessage(&msg);	DispatchMessage(&msg);
	}
	return _close;
}

static int	getkey(void){
	INPUT_RECORD	buf;
	HANDLE	hstdin=GetStdHandle(STD_INPUT_HANDLE);
	do{
		DWORD	n;
		ReadConsoleInput(hstdin,&buf,1,&n);
	}while( buf.Event.KeyEvent.bKeyDown==1);
	FlushConsoleInputBuffer(hstdin);
	return buf.Event.KeyEvent.uChar.AsciiChar;
}

#define	sleep(millsec)	Sleep((DWORD)(millsec))

//////////////	Animation
#define	MEM_MAX	64								// 64:: メモリ配列の最大数
static int	_buffnmb,_mwide,_mhite;					// バファ数。幅 & 高さ
static HDC	_hdc3;
static HBITMAP	_hbitmap[MEM_MAX];

static void	swapbuffers(void){						// ginit(2)にする
	if(_sort==2){
		BitBlt(_hdc2,0,0,_wide,_hite,_hdc,0,0,SRCCOPY);
		gcls();
	}
}

static void	setbitmap(int nm, int wide,int hite){
	int	j;
	if(nm>MEM_MAX){
		printf("フレーム枚数は%d枚までです。", MEM_MAX);	exit(-1);
	}
	if(_sort==3){
		_buffnmb=nm;	_mwide=wide;		_mhite=hite;
		for(j=0;j<nm;j++)
			_hbitmap[j]=CreateCompatibleBitmap(_hdc,_mwide,_mhite);
		_hdc3=CreateCompatibleDC(_hdc);
	}
}

static void	getimage(int n,int x,int y){
	x= x+_xorg;	y= _hite-y-_yorg;
	SelectObject(_hdc3,_hbitmap[n]);
	BitBlt(_hdc3,0,0,_mwide,_mhite,_hdc,x,y,SRCCOPY);
}

static void	putimage(int n,int x,int y){
	x= x+_xorg;	y= _hite-y-_yorg;
	SelectObject(_hdc3,_hbitmap[n]);
	BitBlt(_hdc,x,y,_mwide,_mhite,_hdc3,0,0,SRCCOPY);
}

static void	paintimage(int n,int x,int y){
	x= x+_xorg;	y= _hite-y-_yorg;
	SelectObject(_hdc3,_hbitmap[n]);
	BitBlt(_hdc,x,y,_mwide,_mhite,_hdc3,0,0,SRCPAINT);
}

static void	setlayer(int nm){
	int	j;
	if(nm>MEM_MAX){
		printf("レイヤー数は%d枚までです。", MEM_MAX);	exit(-1);
	}
	if(_sort==2){
		HDC	hdc=CreateCompatibleDC(_hdc2);
		_buffnmb=nm;
		for(j=0;j<nm;j++){
			_hbitmap[j]=CreateCompatibleBitmap(_hdc2,_wide,_hite);
			SelectObject(hdc,_hbitmap[j]);
			PatBlt(hdc,0,0,_wide,_hite,BLACKNESS);
		}
		DeleteDC(hdc);
	}
}

static void	selectlayer(int n){
	SelectObject(_hdc, n<0? _hbmp: _hbitmap[n]);
}

static void	putlayers(void){
	int	j;
	HDC	hdc=CreateCompatibleDC(_hdc2);
	SelectObject(_hdc,_hbmp);
	for(j=0; j<_buffnmb; j++){
		SelectObject(hdc,_hbitmap[j]);
		BitBlt(_hdc,0,0,_wide,_hite,hdc,0,0,SRCPAINT);
	}
	DeleteDC(hdc);
}

///////// Bitmap
static void	_createbitmap(void){
	if(_sort==2){
		_hdc=CreateCompatibleDC(_hdc2);
		_hbmp=CreateCompatibleBitmap(_hdc2,_wide,_hite);
		SelectObject(_hdc,_hbmp);
	}else{
		SetMapMode(_hdc,MM_TEXT);		// 1 or 3
		SelectClipRgn(_hdc,NULL);
		_hdc2=CreateCompatibleDC(_hdc);
		_hbmp=CreateCompatibleBitmap(_hdc,_wide,_hite);
		SelectObject(_hdc2,_hbmp);
	}
}

static void	_delbitmap(HWND hwnd){
	int	j;
	DeleteObject(_hbmp);
	for(j=0;j<_buffnmb;j++)	DeleteObject(_hbitmap[j]);
	if(_hdc3)	DeleteDC(_hdc3);
	if(_sort==2){
		DeleteDC(_hdc);
		ReleaseDC(hwnd, _hdc2);
	}else{
		DeleteDC(_hdc2);						// _sort==1 or 3
		EndPaint(hwnd,&_ghpaint);
	}
}

/////////	Windows API
static LRESULT CALLBACK	_window_proc(HWND hwnd,UINT uMsg,WPARAM wParam,LPARAM lParam){
	static	int	f=0;			// f: flag for WM_CLOSE,no flicker
	switch(uMsg){
		case WM_PAINT:
			ValidateRect(hwnd,NULL);
			if(_sort==2){
				if(f||_end)	BitBlt(_hdc2,0,0,_wide,_hite,_hdc,0,0,SRCCOPY);
			}else	BitBlt(_hdc,0,0,_wide,_hite,_hdc2,0,0,SRCCOPY);
			return 0;
		case WM_KEYDOWN:
		case WM_CLOSE:
			if(_sort==2){
				BitBlt(_hdc,0,0,_wide,_hite,_hdc2,0,0,SRCCOPY);	f=1;
			}else	BitBlt(_hdc2,0,0,_wide,_hite,_hdc,0,0,SRCCOPY);
			if(MessageBox(hwnd,_T("FINISH ?"),_T(""),MB_YESNO|MB_ICONQUESTION)==IDNO){
				f=0;	return 0;
			}
			_close=1;	_delbitmap(hwnd);	DestroyWindow(hwnd);
			return 0;
		case WM_DESTROY:
			PostQuitMessage(0);	return 0;
		default:	return DefWindowProc(hwnd,uMsg,wParam,lParam);
	}
}

static int	_gint(void){
	char	fn[FILENAME_MAX];
	WNDCLASS	wc;
	HWND	hwnd;
	DWORD	wstyle=WS_OVERLAPPEDWINDOW ;
	int	t_h=GetSystemMetrics(SM_CYCAPTION);		// Title height
	int	b_w=GetSystemMetrics(SM_CYFRAME);		// Frame width
	wc.lpszClassName=_T("_gint");
	wc.lpszMenuName=NULL;
	wc.hInstance=GetModuleHandle(NULL);
	wc.lpfnWndProc=_window_proc;		wc.hIcon=NULL;
	wc.hCursor=LoadCursor(NULL,IDC_ARROW);
	wc.hbrBackground=(HBRUSH)GetStockObject(BLACK_BRUSH);
	wc.style=0;		wc.cbClsExtra=0;	wc.cbWndExtra=0;
	if(!RegisterClass(&wc))  return 0;
	GetModuleFileName(NULL,(LPTCH)fn,FILENAME_MAX);
	hwnd=CreateWindow(wc.lpszClassName,(LPTSTR)fn,wstyle,_leftmrg,_topmrg,
		_wide+2*b_w,_hite+t_h+2*b_w,NULL,NULL,wc.hInstance,NULL);
	SetForegroundWindow(hwnd);					// disp on front
	ShowWindow(hwnd,SW_SHOWDEFAULT);
	if(_sort==2){
		_hdc2=GetDC(hwnd);
	}else{
		_hdc=BeginPaint(hwnd,&_ghpaint);				// _sort==1or3
		FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
	}
	return 0;
}

static void	_gend(void){
	MSG	msg;
	if(_sort==2){
		BitBlt(_hdc,0,0,_wide,_hite,_hdc2,0,0,SRCCOPY);	_end=1;
	}else	BitBlt(_hdc2,0,0,_wide,_hite,_hdc,0,0,SRCCOPY);
	while(GetMessage(&msg,NULL,0,0))	DispatchMessage(&msg);
}

static int	ginit(int st){			// st:1 一枚, 2 二枚, 3 アニメ用
	_sort=st;	_gint();	_createbitmap();	atexit(_gend);
	return 0;
}

static void	gwinsize(int w, int h){
	_wide=w;	_hite=h;
}

static void	gwinpos(int left, int top){
	_leftmrg=left;	_topmrg=top;
}

static void	gsetorg(int x, int y){
	_xorg=x;	_yorg=y;
}

/////////	ビットマップファイル入出力	(cm: 1 or 3, datt:one_dim_array)
static int	savebmpfile(char *filnam,int cm,int wd,int ht, unsigned char *datt){
	int	pad=((4-(wd % 4))*cm) % 4,	hdsz=54+(cm==1 ? 1024 : 0);
	long	int	imgsz=(wd*cm+pad)*ht,	fsiz=imgsz+hdsz;
	int	h[54];							// h[54]:BMP_header
	int	x,y,c,k;
	FILE	*fw;
#if defined(_MSC_VER) && _MSC_VER >= 1400		// later VC8.0(VC++2005)
	fopen_s(&fw,filnam,"wb");
#else
	fw=fopen(filnam,"wb");
#endif
	if(fw==NULL){
		printf("Can't Open \"%s\".\n", filnam);	exit(-1);
	}
	for(k=0; k<54; k++)	h[k]=0;							// h[]:initialize
	h[ 0]='B';	h[ 1]='M';
	h[ 2]= fsiz % 256;	h[ 3]=(fsiz>>8) % 256;	h[ 4]=(fsiz>>16) % 256;
	h[ 5]=(fsiz >>24) % 256;		//// h[6], h[7], h[8], h[9]:: reserve
	h[10]= hdsz % 256;	h[11]=(hdsz>>8) % 256;  // h[12] = h[13] = 0
	h[14]= 40;		// h[15], h[16], h[17] :: information size
	h[18]= wd % 256;	h[19]=(wd>>8) % 256;	h[20]=(wd>>16) % 256;
	h[21]=(wd>>24) % 256;
	h[22]= ht % 256;	h[23]=(ht>>8) % 256;	h[24]=(ht>>16) % 256;
	h[25]=(ht>>24) % 256;
	h[26]= 1;										//// pict nmb. h[27]=0
	h[28]= cm*8;				// h[30] ... h[53] use initializing values
	h[34]=imgsz%256;	h[35]=(imgsz>>8)%256;	h[36]=(imgsz>>16)%256;
	h[37]=(imgsz>>24)%256;
	for(k=0; k<54; k++)	fputc(h[k],fw);			// write header & palette
	if(cm==1){
		for(k=0; k<256; k++){						// save palette
			fputc(k,fw);	fputc(k,fw);	fputc(k,fw);	fputc(0,fw);
		}
	}
	for(y=0; y<ht; y++){
		for(x=0; x<wd; x++){							// write pict_data
			for(c=0; c<cm; c++)	fputc(datt[(y*wd+x)*cm+c],fw);
		}
		for(k=0; k<pad; k++) fputc(0,fw);				// write padding
	}
	fclose(fw);
	return 0;
}

static int	loadbmpfile(char *loadfil,int *cm,int *wd,int *ht,unsigned char *dat){
	int	h[54], pal[1024], pad, j, k, c;
	FILE	*fr;
#if defined(_MSC_VER) && _MSC_VER >= 1400		// later VC8.0(VC++2005)
	if(strstr(loadfil,".bmp")==NULL)	strcat_s(loadfil,128,".bmp");
	fopen_s(&fr,loadfil,"rb");
#else
	if(strstr(loadfil,".bmp")==NULL)	strcat(loadfil,".bmp");
	fr=fopen(loadfil,"rb");
#endif
	if(fr==NULL){
		printf("Can't Open \"%s\".\n", loadfil);	exit(-1);
	}
	for(j=0; j<54; j++)	h[j]=fgetc(fr);						// read header
	if(h[0]!='B'&& h[1]!='M'){
		printf("Not a BMP file!\n");	exit(-1);
	}
	*wd=h[18]+(h[19]<<8)+(h[20]<<16)+(h[21]<<24);
	*ht=h[22]+(h[23]<<8)+(h[24]<<16)+(h[25]<<24);
	*cm=h[28]/8;							// cm= 1:mono, 3:full_color
	if(*cm==1)	for(j=0; j<1024; j++)	pal[j]=fgetc(fr);	// read palette
	pad=((4-(*wd % 4))*(*cm)) % 4;					// calc. padding
	for(j=0; j<*ht; j++){							// read pict_data
		for(k=0; k<*wd; k++){
			for(c=0; c<*cm ; c++)	*(dat++)=(unsigned char)fgetc(fr);
		}
		for(k=0; k<pad; k++)	fgetc(fr);			// skip padding
	}
	fclose(fr);
	return 0;
}
///////////	画面の保存
static void	save_screen(char *fname)
{
	int	cm=3, pad=((4-_wide%4)*cm)%4, hdsz=54+(cm==1?1024:0);
	long	int	imgsz=(_wide*cm+pad)*_hite, fsiz=imgsz+hdsz;
	BITMAPFILEHEADER	bmf={'B'+('M'<<8)};
	BITMAPINFO		bmi={40};
	unsigned	char	*data=(unsigned char*)malloc(imgsz);
	FILE	*fp;
	HDC		hdc =CreateCompatibleDC(_hdc);
	HBITMAP	hbmp=CreateCompatibleBitmap(_hdc,_wide,_hite),
			hbmp_old=(HBITMAP)SelectObject(hdc,hbmp);
	bmf.bfSize=fsiz;
	bmf.bfOffBits=hdsz;
	bmi.bmiHeader.biWidth=_wide;
	bmi.bmiHeader.biHeight=_hite;
	bmi.bmiHeader.biPlanes=1;
	bmi.bmiHeader.biBitCount=(short)(cm*8);
	bmi.bmiHeader.biSizeImage=imgsz;
	BitBlt(hdc,0,0,_wide,_hite,_hdc,0,0,SRCCOPY);
	GetDIBits(_hdc, hbmp, 0, _hite, data, &bmi, DIB_RGB_COLORS);
	SelectObject(hdc, hbmp_old);	DeleteObject(hbmp);	DeleteDC(hdc);
#if defined(_MSC_VER) && _MSC_VER >= 1400	// VC8.0(VC++2005)以降なら
	fopen_s(&fp,fname,"wb");
#else
	fp=fopen(fname,"wb");
#endif
	fwrite(&bmf,14,1,fp);
	fwrite(&bmi.bmiHeader, bmi.bmiHeader.biSize, 1, fp);
	fwrite(data,imgsz,1,fp);
	fclose(fp);
	free(data);
	printf("画像を\"%s\"に保存しました。\n", fname);
}

#endif // ndef _WINGXA_H
