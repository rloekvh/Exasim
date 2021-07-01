template <typename T> void opuUbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T odg1 = odg[0*ng+i];
		T odg2 = odg[1*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = 1.0/udg1;
		T t3 = 1.0/(udg1*udg1);
		T t4 = udg2*udg2;
		T t5 = t3*t4;
		T t6 = udg3*udg3;
		T t7 = t3*t6;
		T t8 = t5+t7;
		T t10 = (t8*udg1)/2.0;
		T t9 = -t10+udg4;
		T t11 = param1-1.0;
		T t12 = nlg1*udg2;
		T t13 = nlg2*udg3;
		T t14 = t12+t13;
		T t15 = param1*t2*t9*t11;
		T t16 = sqrt(t15);
		T t17 = 1.0/param1;
		T t18 = t2*t14*1.0E3;
		T t19 = 1.0/t9;
		T t20 = t16*1.0E3;
		T t21 = t2*udg4;
		T t22 = t2*t9*t11;
		T t23 = t2*t14;
		T t24 = 1.0/t11;
		T t25 = tanh(t18);
		T t26 = t18-t20;
		T t27 = tanh(t26);
		T t28 = nlg1*udg3;
		T t32 = nlg2*udg2;
		T t29 = t28-t32;
		T t30 = t18+t20;
		T t31 = tanh(t30);
		T t33 = t16*t24;
		T t34 = (t17*t19*t27*t29)/2.0;
		T t35 = (t17*t19*t29*t31)/2.0;
		T t45 = t17*t19*t25*t29;
		T t36 = t34+t35-t45;
		T t37 = t23+t33;
		T t38 = (t17*t19*t27*t37*udg1)/2.0;
		T t39 = t23-t33;
		T t40 = (t17*t19*t31*t39*udg1)/2.0;
		T t42 = t14*t17*t19*t25;
		T t41 = t38+t40-t42;
		T t43 = t16-t23;
		T t44 = t16+t23;
		T t46 = nlg2*t36;
		T t47 = param8-udg4;
		T t48 = t29*t29;
		T t49 = t2*t14*t17*t19*t25*t29;
		T t50 = param6-udg2;
		T t51 = t2*t17*t19*t25*t48;
		T t79 = (t2*t17*t19*t27*t48)/2.0;
		T t80 = (t2*t17*t19*t31*t48)/2.0;
		T t52 = t25+t51-t79-t80;
		T t53 = nlg2*t52;
		T t54 = (t17*t19*t27*t29*t43)/2.0;
		T t81 = (t17*t19*t29*t31*t44)/2.0;
		T t55 = t49+t54-t81;
		T t56 = t53-nlg1*t55;
		T t57 = t14*t14;
		T t58 = t2*t17*t19*t25*t57;
		T t59 = (t17*t19*t27*t37*t43*udg1)/2.0;
		T t82 = (t17*t19*t31*t39*t44*udg1)/2.0;
		T t60 = t58+t59-t82;
		T t61 = nlg1*t60;
		T t62 = (t17*t19*t27*t29*t37)/2.0;
		T t63 = (t17*t19*t29*t31*t39)/2.0;
		T t64 = -t49+t62+t63;
		T t65 = nlg2*t64;
		T t66 = t61+t65;
		T t67 = param7-udg3;
		T t68 = param5-udg1;
		T t69 = t2*udg4*2.0;
		T t70 = t2*t9*t11*2.0;
		T t74 = param1*t2*t9*4.0;
		T t71 = t69+t70-t74;
		T t72 = -t16+t23;
		T t77 = t16*t24*t44;
		T t73 = t21+t22-t77;
		T t75 = t16*t24*t72;
		T t76 = t21+t22+t75;
		T t78 = nlg1*t36;
		T t83 = (t17*t19*t27*t29*t72)/2.0;
		T t84 = -t49+t81+t83;
		T t85 = nlg1*t52-nlg2*t84;
		T t86 = (t17*t19*t27*t37*t72*udg1)/2.0;
		T t87 = -t58+t82+t86;
		T t88 = nlg2*t87;
		T t89 = nlg1*t64;
		T t90 = t88+t89;
		T t91 = (t17*t19*t27*t72*t76*udg1)/2.0;
		T t92 = (t17*t19*t31*t44*t73*udg1)/2.0;
		T t93 = t91+t92-(t14*t17*t19*t25*t71)/2.0;
		T t94 = t2*t25*t29;
		T t95 = (t17*t19*t25*t29*t71)/2.0;
		T t96 = t94+t95-(t17*t19*t27*t29*t76)/2.0-(t17*t19*t29*t31*t73)/2.0;
		T t97 = t2*t14*t16;
		T t98 = (t3*t57)/2.0;
		T t99 = (t3*t48)/2.0;
		T t100 = t98+t99;
		T t101 = t21+t22-t97;
		T t102 = t21+t22+t97;
		T t103 = t17*t19*t25*t29*t100;
		T t104 = t94+t103-(t17*t19*t27*t29*t101)/2.0-(t17*t19*t29*t31*t102)/2.0;
		T t105 = (t17*t19*t27*t37*t101*udg1)/2.0;
		T t106 = (t17*t19*t31*t39*t102*udg1)/2.0;
		T t107 = t105+t106-t14*t17*t19*t25*t100;
		f[0*ng+i] = param5/2.0+udg1/2.0-(t68*((t17*t19*t27*udg1*(t21+t22-t16*t24*t43))/2.0-(t17*t19*t25*t71*udg1)/2.0+(t17*t19*t31*t73*udg1)/2.0))/2.0-(t47*(-t17*t19*t25*udg1+(t17*t19*t27*udg1)/2.0+(t17*t19*t31*udg1)/2.0))/2.0-(t50*(t46-nlg1*t41))/2.0+(t67*(t78+nlg2*t41))/2.0;
		f[1*ng+i] = param6/2.0+udg2/2.0+(t47*(t46+nlg1*(t42+(t17*t19*t27*t43*udg1)/2.0-(t17*t19*t31*t44*udg1)/2.0)))/2.0-(t50*(nlg2*t56+nlg1*t66))/2.0+(t67*(nlg1*t56-nlg2*t66))/2.0-(t68*(nlg1*t93+nlg2*t96))/2.0;
		f[2*ng+i] = param7/2.0+udg3/2.0+(t50*(nlg2*t85+nlg1*t90))/2.0-(t67*(nlg1*t85-nlg2*t90))/2.0-(t68*(nlg2*t93-nlg1*t96))/2.0-(t47*(t78+nlg2*(-t42+(t17*t19*t31*t44*udg1)/2.0+(t17*t19*t27*t72*udg1)/2.0)))/2.0;
		f[3*ng+i] = param8/2.0+udg4/2.0-(t47*(-t17*t19*t25*t100*udg1+(t17*t19*t27*t101*udg1)/2.0+(t17*t19*t31*t102*udg1)/2.0))/2.0+(t68*(t3*t25*t48+(t17*t19*t25*t71*t100*udg1)/2.0-(t17*t19*t27*t76*t101*udg1)/2.0-(t17*t19*t31*t73*t102*udg1)/2.0))/2.0+(t50*(nlg2*t104+nlg1*t107))/2.0-(t67*(nlg1*t104-nlg2*t107))/2.0;
	}
}

template <typename T> void opuUbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T odg1 = odg[0*ng+i];
		T odg2 = odg[1*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = udg4;
	}
}

template <typename T> void opuUbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T odg1 = odg[0*ng+i];
		T odg2 = odg[1*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = (param9*param11*udg1)/param10;
	}
}

template <typename T> void opuUbou4(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T odg1 = odg[0*ng+i];
		T odg2 = odg[1*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg2;
		T t3 = nlg2*udg3;
		T t4 = t2+t3;
		f[0*ng+i] = udg1;
		f[1*ng+i] = udg2-nlg1*t4;
		f[2*ng+i] = udg3-nlg2*t4;
		f[3*ng+i] = udg4;
	}
}

template <typename T> void opuUbou5(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T odg1 = odg[0*ng+i];
		T odg2 = odg[1*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		f[0*ng+i] = param5;
		f[1*ng+i] = param6;
		f[2*ng+i] = param7;
		f[3*ng+i] = param8;
	}
}

template <typename T> void opuUbou6(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T odg1 = odg[0*ng+i];
		T odg2 = odg[1*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = udg2;
		f[2*ng+i] = udg3;
		f[3*ng+i] = udg4;
	}
}

template <typename T> void opuUbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (ib == 1)
		opuUbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		opuUbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		opuUbou3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		opuUbou4(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		opuUbou5(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 6)
		opuUbou6(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void opuUbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void opuUbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
