/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6058223124271342433);
void inv_err_fun(double *nom_x, double *true_x, double *out_4594994462099967109);
void H_mod_fun(double *state, double *out_7455493487645305107);
void f_fun(double *state, double dt, double *out_3779774445990277428);
void F_fun(double *state, double dt, double *out_9123363132626355535);
void h_3(double *state, double *unused, double *out_4921602961994071490);
void H_3(double *state, double *unused, double *out_7842695018323111031);
void h_4(double *state, double *unused, double *out_863731571442044149);
void H_4(double *state, double *unused, double *out_7157404151381093669);
void h_9(double *state, double *unused, double *out_421674556503941375);
void H_9(double *state, double *unused, double *out_4834898607510237699);
void h_10(double *state, double *unused, double *out_9210518317766081543);
void H_10(double *state, double *unused, double *out_3566185315675630910);
void h_12(double *state, double *unused, double *out_7575109034293337248);
void H_12(double *state, double *unused, double *out_5050848639641884100);
void h_31(double *state, double *unused, double *out_3085913874242626197);
void H_31(double *state, double *unused, double *out_7725938414785303921);
void h_32(double *state, double *unused, double *out_2771588110458732063);
void H_32(double *state, double *unused, double *out_8770303937276496048);
void h_13(double *state, double *unused, double *out_1115330517262319623);
void H_13(double *state, double *unused, double *out_4526753010816576093);
void h_14(double *state, double *unused, double *out_421674556503941375);
void H_14(double *state, double *unused, double *out_4834898607510237699);
void h_19(double *state, double *unused, double *out_2905258556630152314);
void H_19(double *state, double *unused, double *out_6643810532227050307);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);