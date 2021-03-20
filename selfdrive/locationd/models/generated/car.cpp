
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2622791219461976747) {
   out_2622791219461976747[0] = delta_x[0] + nom_x[0];
   out_2622791219461976747[1] = delta_x[1] + nom_x[1];
   out_2622791219461976747[2] = delta_x[2] + nom_x[2];
   out_2622791219461976747[3] = delta_x[3] + nom_x[3];
   out_2622791219461976747[4] = delta_x[4] + nom_x[4];
   out_2622791219461976747[5] = delta_x[5] + nom_x[5];
   out_2622791219461976747[6] = delta_x[6] + nom_x[6];
   out_2622791219461976747[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7195422080887650603) {
   out_7195422080887650603[0] = -nom_x[0] + true_x[0];
   out_7195422080887650603[1] = -nom_x[1] + true_x[1];
   out_7195422080887650603[2] = -nom_x[2] + true_x[2];
   out_7195422080887650603[3] = -nom_x[3] + true_x[3];
   out_7195422080887650603[4] = -nom_x[4] + true_x[4];
   out_7195422080887650603[5] = -nom_x[5] + true_x[5];
   out_7195422080887650603[6] = -nom_x[6] + true_x[6];
   out_7195422080887650603[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_1888720099587935626) {
   out_1888720099587935626[0] = 1.0;
   out_1888720099587935626[1] = 0.0;
   out_1888720099587935626[2] = 0.0;
   out_1888720099587935626[3] = 0.0;
   out_1888720099587935626[4] = 0.0;
   out_1888720099587935626[5] = 0.0;
   out_1888720099587935626[6] = 0.0;
   out_1888720099587935626[7] = 0.0;
   out_1888720099587935626[8] = 0.0;
   out_1888720099587935626[9] = 1.0;
   out_1888720099587935626[10] = 0.0;
   out_1888720099587935626[11] = 0.0;
   out_1888720099587935626[12] = 0.0;
   out_1888720099587935626[13] = 0.0;
   out_1888720099587935626[14] = 0.0;
   out_1888720099587935626[15] = 0.0;
   out_1888720099587935626[16] = 0.0;
   out_1888720099587935626[17] = 0.0;
   out_1888720099587935626[18] = 1.0;
   out_1888720099587935626[19] = 0.0;
   out_1888720099587935626[20] = 0.0;
   out_1888720099587935626[21] = 0.0;
   out_1888720099587935626[22] = 0.0;
   out_1888720099587935626[23] = 0.0;
   out_1888720099587935626[24] = 0.0;
   out_1888720099587935626[25] = 0.0;
   out_1888720099587935626[26] = 0.0;
   out_1888720099587935626[27] = 1.0;
   out_1888720099587935626[28] = 0.0;
   out_1888720099587935626[29] = 0.0;
   out_1888720099587935626[30] = 0.0;
   out_1888720099587935626[31] = 0.0;
   out_1888720099587935626[32] = 0.0;
   out_1888720099587935626[33] = 0.0;
   out_1888720099587935626[34] = 0.0;
   out_1888720099587935626[35] = 0.0;
   out_1888720099587935626[36] = 1.0;
   out_1888720099587935626[37] = 0.0;
   out_1888720099587935626[38] = 0.0;
   out_1888720099587935626[39] = 0.0;
   out_1888720099587935626[40] = 0.0;
   out_1888720099587935626[41] = 0.0;
   out_1888720099587935626[42] = 0.0;
   out_1888720099587935626[43] = 0.0;
   out_1888720099587935626[44] = 0.0;
   out_1888720099587935626[45] = 1.0;
   out_1888720099587935626[46] = 0.0;
   out_1888720099587935626[47] = 0.0;
   out_1888720099587935626[48] = 0.0;
   out_1888720099587935626[49] = 0.0;
   out_1888720099587935626[50] = 0.0;
   out_1888720099587935626[51] = 0.0;
   out_1888720099587935626[52] = 0.0;
   out_1888720099587935626[53] = 0.0;
   out_1888720099587935626[54] = 1.0;
   out_1888720099587935626[55] = 0.0;
   out_1888720099587935626[56] = 0.0;
   out_1888720099587935626[57] = 0.0;
   out_1888720099587935626[58] = 0.0;
   out_1888720099587935626[59] = 0.0;
   out_1888720099587935626[60] = 0.0;
   out_1888720099587935626[61] = 0.0;
   out_1888720099587935626[62] = 0.0;
   out_1888720099587935626[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_9052397247224008513) {
   out_9052397247224008513[0] = state[0];
   out_9052397247224008513[1] = state[1];
   out_9052397247224008513[2] = state[2];
   out_9052397247224008513[3] = state[3];
   out_9052397247224008513[4] = state[4];
   out_9052397247224008513[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_9052397247224008513[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_9052397247224008513[7] = state[7];
}
void F_fun(double *state, double dt, double *out_5337452365801478129) {
   out_5337452365801478129[0] = 1;
   out_5337452365801478129[1] = 0;
   out_5337452365801478129[2] = 0;
   out_5337452365801478129[3] = 0;
   out_5337452365801478129[4] = 0;
   out_5337452365801478129[5] = 0;
   out_5337452365801478129[6] = 0;
   out_5337452365801478129[7] = 0;
   out_5337452365801478129[8] = 0;
   out_5337452365801478129[9] = 1;
   out_5337452365801478129[10] = 0;
   out_5337452365801478129[11] = 0;
   out_5337452365801478129[12] = 0;
   out_5337452365801478129[13] = 0;
   out_5337452365801478129[14] = 0;
   out_5337452365801478129[15] = 0;
   out_5337452365801478129[16] = 0;
   out_5337452365801478129[17] = 0;
   out_5337452365801478129[18] = 1;
   out_5337452365801478129[19] = 0;
   out_5337452365801478129[20] = 0;
   out_5337452365801478129[21] = 0;
   out_5337452365801478129[22] = 0;
   out_5337452365801478129[23] = 0;
   out_5337452365801478129[24] = 0;
   out_5337452365801478129[25] = 0;
   out_5337452365801478129[26] = 0;
   out_5337452365801478129[27] = 1;
   out_5337452365801478129[28] = 0;
   out_5337452365801478129[29] = 0;
   out_5337452365801478129[30] = 0;
   out_5337452365801478129[31] = 0;
   out_5337452365801478129[32] = 0;
   out_5337452365801478129[33] = 0;
   out_5337452365801478129[34] = 0;
   out_5337452365801478129[35] = 0;
   out_5337452365801478129[36] = 1;
   out_5337452365801478129[37] = 0;
   out_5337452365801478129[38] = 0;
   out_5337452365801478129[39] = 0;
   out_5337452365801478129[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5337452365801478129[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5337452365801478129[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5337452365801478129[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5337452365801478129[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5337452365801478129[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5337452365801478129[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5337452365801478129[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5337452365801478129[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5337452365801478129[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5337452365801478129[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5337452365801478129[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5337452365801478129[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5337452365801478129[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5337452365801478129[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5337452365801478129[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5337452365801478129[56] = 0;
   out_5337452365801478129[57] = 0;
   out_5337452365801478129[58] = 0;
   out_5337452365801478129[59] = 0;
   out_5337452365801478129[60] = 0;
   out_5337452365801478129[61] = 0;
   out_5337452365801478129[62] = 0;
   out_5337452365801478129[63] = 1;
}
void h_25(double *state, double *unused, double *out_3216024780874029830) {
   out_3216024780874029830[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7806941210924074381) {
   out_7806941210924074381[0] = 0;
   out_7806941210924074381[1] = 0;
   out_7806941210924074381[2] = 0;
   out_7806941210924074381[3] = 0;
   out_7806941210924074381[4] = 0;
   out_7806941210924074381[5] = 0;
   out_7806941210924074381[6] = 1;
   out_7806941210924074381[7] = 0;
}
void h_24(double *state, double *unused, double *out_4742376711754663176) {
   out_4742376711754663176[0] = state[4];
   out_4742376711754663176[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3707274211322865014) {
   out_3707274211322865014[0] = 0;
   out_3707274211322865014[1] = 0;
   out_3707274211322865014[2] = 0;
   out_3707274211322865014[3] = 0;
   out_3707274211322865014[4] = 1;
   out_3707274211322865014[5] = 0;
   out_3707274211322865014[6] = 0;
   out_3707274211322865014[7] = 0;
   out_3707274211322865014[8] = 0;
   out_3707274211322865014[9] = 0;
   out_3707274211322865014[10] = 0;
   out_3707274211322865014[11] = 0;
   out_3707274211322865014[12] = 0;
   out_3707274211322865014[13] = 1;
   out_3707274211322865014[14] = 0;
   out_3707274211322865014[15] = 0;
}
void h_30(double *state, double *unused, double *out_3554810445476321106) {
   out_3554810445476321106[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8852986988659965724) {
   out_8852986988659965724[0] = 0;
   out_8852986988659965724[1] = 0;
   out_8852986988659965724[2] = 0;
   out_8852986988659965724[3] = 0;
   out_8852986988659965724[4] = 1;
   out_8852986988659965724[5] = 0;
   out_8852986988659965724[6] = 0;
   out_8852986988659965724[7] = 0;
}
void h_26(double *state, double *unused, double *out_3447952258316744673) {
   out_3447952258316744673[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6465505077228426829) {
   out_6465505077228426829[0] = 0;
   out_6465505077228426829[1] = 0;
   out_6465505077228426829[2] = 0;
   out_6465505077228426829[3] = 0;
   out_6465505077228426829[4] = 0;
   out_6465505077228426829[5] = 0;
   out_6465505077228426829[6] = 0;
   out_6465505077228426829[7] = 1;
}
void h_27(double *state, double *unused, double *out_2272514927972367395) {
   out_2272514927972367395[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8306175097212960580) {
   out_8306175097212960580[0] = 0;
   out_8306175097212960580[1] = 0;
   out_8306175097212960580[2] = 0;
   out_8306175097212960580[3] = 1;
   out_8306175097212960580[4] = 0;
   out_8306175097212960580[5] = 0;
   out_8306175097212960580[6] = 0;
   out_8306175097212960580[7] = 0;
}
void h_29(double *state, double *unused, double *out_2691640611983852362) {
   out_2691640611983852362[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1232399086169572724) {
   out_1232399086169572724[0] = 0;
   out_1232399086169572724[1] = 1;
   out_1232399086169572724[2] = 0;
   out_1232399086169572724[3] = 0;
   out_1232399086169572724[4] = 0;
   out_1232399086169572724[5] = 0;
   out_1232399086169572724[6] = 0;
   out_1232399086169572724[7] = 0;
}
void h_28(double *state, double *unused, double *out_5814375002743371124) {
   out_5814375002743371124[0] = state[5];
   out_5814375002743371124[1] = state[6];
}
void H_28(double *state, double *unused, double *out_7708915649897997980) {
   out_7708915649897997980[0] = 0;
   out_7708915649897997980[1] = 0;
   out_7708915649897997980[2] = 0;
   out_7708915649897997980[3] = 0;
   out_7708915649897997980[4] = 0;
   out_7708915649897997980[5] = 1;
   out_7708915649897997980[6] = 0;
   out_7708915649897997980[7] = 0;
   out_7708915649897997980[8] = 0;
   out_7708915649897997980[9] = 0;
   out_7708915649897997980[10] = 0;
   out_7708915649897997980[11] = 0;
   out_7708915649897997980[12] = 0;
   out_7708915649897997980[13] = 0;
   out_7708915649897997980[14] = 1;
   out_7708915649897997980[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
