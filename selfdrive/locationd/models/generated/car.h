/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2622791219461976747);
void inv_err_fun(double *nom_x, double *true_x, double *out_7195422080887650603);
void H_mod_fun(double *state, double *out_1888720099587935626);
void f_fun(double *state, double dt, double *out_9052397247224008513);
void F_fun(double *state, double dt, double *out_5337452365801478129);
void h_25(double *state, double *unused, double *out_3216024780874029830);
void H_25(double *state, double *unused, double *out_7806941210924074381);
void h_24(double *state, double *unused, double *out_4742376711754663176);
void H_24(double *state, double *unused, double *out_3707274211322865014);
void h_30(double *state, double *unused, double *out_3554810445476321106);
void H_30(double *state, double *unused, double *out_8852986988659965724);
void h_26(double *state, double *unused, double *out_3447952258316744673);
void H_26(double *state, double *unused, double *out_6465505077228426829);
void h_27(double *state, double *unused, double *out_2272514927972367395);
void H_27(double *state, double *unused, double *out_8306175097212960580);
void h_29(double *state, double *unused, double *out_2691640611983852362);
void H_29(double *state, double *unused, double *out_1232399086169572724);
void h_28(double *state, double *unused, double *out_5814375002743371124);
void H_28(double *state, double *unused, double *out_7708915649897997980);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
