////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005-2022 The Octave Project Developers
// Copyright (C) 2022 Julian Hall
//
// See the file COPYRIGHT.md in the top-level directory of this
// distribution or <https://octave.org/copyright/>.
//
// This file is part of Octave.
//
// Octave is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, see
// <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////

#include <ctime>

#include <limits>

#include <iostream>

#include <octave/oct.h>

#include <Highs.h>

struct control_params
{
  int msglev;
  int dual;
  int price;
  int itlim;
  int outfrq;
  int branch;
  int btrack;
  int presol;
  int rtest;
  int tmlim;
  int outdly;
  double tolbnd;
  double toldj;
  double tolpiv;
  double objll;
  double objul;
  double tolint;
  double tolobj;
};

const bool dev_run = false;
// GLPK variable types
const int GLP_CV = 1;
const int GLP_IV = 2;
const int GLP_BV = 3;

static int
highs (int sense, int n, int m, double *c, int nz, int *rn, int *cn,
      double *a, double *b, char *ctype, int *freeLB, double *lb,
      int *freeUB, double *ub, int *vartype, int isMIP, int lpsolver,
      int save_pb, int scale, const control_params& par,
      double *xmin, double& fmin, int& status,
      double *lambda, double *redcosts, double& time)
{
  std::cout << "Running HiGHS\n";

  int errnum = 0;

  time = 0.0;
  status = -1;    // Initialize status to "bad" value

  clock_t t_start = clock ();

  HighsLp lp;
  const double inf = kHighsInf;

  if (sense != 1) lp.sense_ = ObjSense::kMaximize;

  lp.num_col_ = n;
  lp.num_row_ = m;
  for (int i = 0; i < n; i++) {
    lp.col_cost_.push_back(c[i]);
    lp.col_lower_.push_back(lb[i]);
    lp.col_upper_.push_back(ub[i]);
    if (isMIP) {
      if (vartype[i] == GLP_CV) {
        lp.integrality_.push_back(HighsVarType::kContinuous);
      } else {
        lp.integrality_.push_back(HighsVarType::kInteger);
      }
    }
  }
  for (int i = 0; i < m; i++) {
    double lower = -inf;
    double upper = inf;
    switch (ctype[i])
      {
      case 'F':
        break;
        
      case 'U':
        upper = b[i];
        break;
        
      case 'L':
        lower = b[i];
        break;
        
      case 'S':
        upper = b[i];
        lower = b[i];
        break;
        
      case 'D':
        upper = b[i];
        lower = b[i];
        assert("Row is boxed" == "");
        break;
      }
    lp.row_lower_.push_back(lower);
    lp.row_upper_.push_back(upper);
  }

  std::cout << "Defined consts and bounds\n";
  // Convert triplet form of constraint matrix into CSC
  std::vector<int>& start = lp.a_matrix_.start_;
  std::vector<int> length;
  length.assign(m, 0);
  lp.a_matrix_.index_.resize(nz);
  lp.a_matrix_.value_.resize(nz);
  for (int k = 1; k <= nz; k++) length[cn[k]-1]++;
  for (int i = 0; i < n; i++) {
    start.push_back(start[i] + length[i]);
    length[i] = start[i];
  }
  for (int k = 1; k <= nz; k++) {
    int j = cn[k]-1;
    lp.a_matrix_.index_[length[j]] = rn[k]-1;
    lp.a_matrix_.value_[length[j]] = a[k];
    length[j]++;
  }
  assert(length[n-1] == nz);
  std::cout << "Defined LP\n";
  Highs highs;
  HighsStatus return_status;
  if (!dev_run && par.msglev == 0) {
    return_status = highs.setOptionValue("output_flag", false);
    assert(return_status == HighsStatus::kOk);
  }
  if (par.dual == 0) {
    return_status = highs.setOptionValue("simplex_strategy", 4);
    assert(return_status == HighsStatus::kOk);
  }
  if (par.price == 17) {
    return_status = highs.setOptionValue("simplex_dual_edge_weight_strategy", 0);
    assert(return_status == HighsStatus::kOk);
    return_status = highs.setOptionValue("simplex_primal_edge_weight_strategy", 0);
    assert(return_status == HighsStatus::kOk);
  }
  if (par.itlim < kHighsIInf) {
    return_status = highs.setOptionValue("simplex_iteration_limit", par.itlim);
    assert(return_status == HighsStatus::kOk);
  }
  if (par.presol != 1) {
    return_status = highs.setOptionValue("presolve", kHighsOffString);
    assert(return_status == HighsStatus::kOk);
  }
  if (par.tmlim < kHighsInf) {
    return_status = highs.setOptionValue("time_limit", (double)par.tmlim);
    assert(return_status == HighsStatus::kOk);
  }
  return_status = highs.setOptionValue("primal_feasibility_tolerance", par.tolbnd);
  assert(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("dual_feasibility_tolerance", par.toldj);
  assert(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("mip_feasibility_tolerance", par.tolint);
  assert(return_status == HighsStatus::kOk);
  if (lpsolver == 2) {
    return_status = highs.setOptionValue("solver", kIpmString);
    assert(return_status == HighsStatus::kOk);
  }

  return_status = highs.passModel(lp);
  std::cout << "Pass model yields "<< (int)return_status << "\n";
  assert(return_status == HighsStatus::kOk);

  if (save_pb) {
    return_status = highs.writeModel("outpb.lp");
    assert(return_status == HighsStatus::kOk);
  }

  return_status = highs.run();
  std::cout << "Run model yields "<< (int)return_status << "\n";
  assert(return_status == HighsStatus::kOk);
  const HighsInfo& info = highs.getInfo();
  fmin = info.objective_function_value;
  const HighsSolution& solution = highs.getSolution();
  if (solution.value_valid) {
    for (int i = 0; i < n; i++)
      xmin[i] = solution.col_value[i];
  }
  if (solution.dual_valid) {
    for (int i = 0; i < n; i++)
      redcosts[i] = solution.col_dual[i];
    for (int i = 0; i < m; i++)
      lambda[i] = solution.row_dual[i];
  }
  HighsModelStatus model_status = highs.getModelStatus();
  
  status = 1;
  if (model_status == HighsModelStatus::kOptimal) {
    status = 5;
  } else if (model_status == HighsModelStatus::kInfeasible) {
    status = 3;
  } else if (model_status == HighsModelStatus::kUnbounded) {
    status = 6;
  } else if (info.primal_solution_status == kSolutionStatusInfeasible) {
    status = 2;
  }
  errnum = 0;
  if (model_status == HighsModelStatus::kSolveError) {
    errnum = 5;
  } if (model_status == HighsModelStatus::kIterationLimit) {
    errnum = 8;
  } if (model_status == HighsModelStatus::kTimeLimit) {
    errnum = 9;
  } if (info.mip_gap > 0) {
    errnum = 14;
  }

  time = (clock () - t_start) / CLOCKS_PER_SEC;

  return errnum;
}

OCTAVE_NAMESPACE_BEGIN

#define OCTAVE_HIGHS_GET_REAL_PARAM(NAME, VAL)                           \
  do                                                                    \
    {                                                                   \
      octave_value tmp = PARAM.getfield (NAME);                         \
                                                                        \
      if (tmp.is_defined ())                                            \
        {                                                               \
          if (! tmp.isempty ())                                        \
            VAL = tmp.xscalar_value ("highs: invalid value in PARAM" NAME); \
          else                                                          \
            error ("highs: invalid value in PARAM" NAME);                \
        }                                                               \
    }                                                                   \
  while (0)

#define OCTAVE_HIGHS_GET_INT_PARAM(NAME, VAL)                            \
  do                                                                    \
    {                                                                   \
      octave_value tmp = PARAM.getfield (NAME);                         \
                                                                        \
      if (tmp.is_defined ())                                            \
        {                                                               \
          if (! tmp.isempty ())                                        \
            VAL = tmp.xint_value ("highs: invalid value in PARAM" NAME); \
          else                                                          \
            error ("highs: invalid value in PARAM" NAME);                \
        }                                                               \
    }                                                                   \
  while (0)

DEFUN_DLD (__highs__, args, ,
           "-*- texinfo -*-\n\
@deftypefn {} {[@var{values}] =} __highs__ (@var{args})\n\
Undocumented internal function.\n\
@end deftypefn")
{
  // FIXME: Should we even need checking for an internal function?
  if (args.length () != 9)
    print_usage ();

  // 1st Input.  A column array containing the objective function coefficients.
  int mrowsc = args(0).rows ();

  Matrix C = args(0).xmatrix_value ("__highs__: invalid value of C");

  double *c = C.fortran_vec ();
  Array<int> rn;
  Array<int> cn;
  ColumnVector a;
  int mrowsA;
  int nz = 0;

  // 2nd Input.  A matrix containing the constraints coefficients.
  // If matrix A is NOT a sparse matrix
  if (args(1).issparse ())
    {
      SparseMatrix A = args(1).xsparse_matrix_value ("__highs__: invalid value of A");

      mrowsA = A.rows ();
      octave_idx_type Anc = A.cols ();
      octave_idx_type Anz = A.nnz ();
      rn.resize (dim_vector (Anz+1, 1));
      cn.resize (dim_vector (Anz+1, 1));
      a.resize (Anz+1, 0.0);

      if (Anc != mrowsc)
        error ("__highs__: invalid value of A");

      for (octave_idx_type j = 0; j < Anc; j++)
        for (octave_idx_type i = A.cidx (j); i < A.cidx (j+1); i++)
          {
            nz++;
            rn(nz) = A.ridx (i) + 1;
            cn(nz) = j + 1;
            a(nz) = A.data(i);
          }
    }
  else
    {
      Matrix A = args(1).xmatrix_value ("__highs__: invalid value of A");

      mrowsA = A.rows ();
      rn.resize (dim_vector (mrowsA*mrowsc+1, 1));
      cn.resize (dim_vector (mrowsA*mrowsc+1, 1));
      a.resize (mrowsA*mrowsc+1, 0.0);

      for (int i = 0; i < mrowsA; i++)
        {
          for (int j = 0; j < mrowsc; j++)
            {
              if (A(i, j) != 0)
                {
                  nz++;
                  rn(nz) = i + 1;
                  cn(nz) = j + 1;
                  a(nz) = A(i, j);
                }
            }
        }

    }

  // 3rd Input.  A column array containing the right-hand side value
  //             for each constraint in the constraint matrix.
  Matrix B = args(2).xmatrix_value ("__highs__: invalid value of B");

  double *b = B.fortran_vec ();

  // 4th Input.  An array of length mrowsc containing the lower
  //             bound on each of the variables.
  Matrix LB = args(3).xmatrix_value ("__highs__: invalid value of LB");

  if (LB.numel () < mrowsc)
    error ("__highs__: invalid dimensions for LB");

  double *lb = LB.fortran_vec ();

  // LB argument, default: Free
  Array<int> freeLB (dim_vector (mrowsc, 1));
  for (int i = 0; i < mrowsc; i++)
    {
      if (math::isinf (lb[i]))
        {
          freeLB(i) = 1;
          lb[i] = -numeric_limits<double>::Inf ();
        }
      else
        freeLB(i) = 0;
    }

  // 5th Input.  An array of at least length numcols containing the upper
  //             bound on each of the variables.
  Matrix UB = args(4).xmatrix_value ("__highs__: invalid value of UB");

  if (UB.numel () < mrowsc)
    error ("__highs__: invalid dimensions for UB");

  double *ub = UB.fortran_vec ();

  Array<int> freeUB (dim_vector (mrowsc, 1));
  for (int i = 0; i < mrowsc; i++)
    {
      if (math::isinf (ub[i]))
        {
          freeUB(i) = 1;
          ub[i] = numeric_limits<double>::Inf ();
        }
      else
        freeUB(i) = 0;
    }

  // 6th Input.  A column array containing the sense of each constraint
  //             in the constraint matrix.
  charMatrix CTYPE = args(5).xchar_matrix_value ("__highs__: invalid value of CTYPE");

  char *ctype = CTYPE.fortran_vec ();

  // 7th Input.  A column array containing the types of the variables.
  charMatrix VTYPE = args(6).xchar_matrix_value ("__highs__: invalid value of VARTYPE");

  Array<int> vartype (dim_vector (mrowsc, 1));
  int isMIP = 0;
  for (int i = 0; i < mrowsc ; i++)
    {
      if (VTYPE(i, 0) == 'I')
        {
          isMIP = 1;
          vartype(i) = GLP_IV;
        }
      else
        vartype(i) = GLP_CV;
    }

  // 8th Input.  Sense of optimization.
  int sense;
  double SENSE = args(7).xscalar_value ("__highs__: invalid value of SENSE");

  if (SENSE >= 0)
    sense = 1;
  else
    sense = -1;

  // 9th Input.  A structure containing the control parameters.
  octave_scalar_map PARAM = args(8).xscalar_map_value ("__highs__: invalid value of PARAM");

  control_params par;

  // Integer parameters

  // Level of messages output by the solver
  par.msglev = 1;
  OCTAVE_HIGHS_GET_INT_PARAM ("msglev", par.msglev);
  if (par.msglev < 0 || par.msglev > 3)
    error ("__highs__: PARAM.msglev must be 0 (no output) or 1 (error and warning messages only [default]) or 2 (normal output) or 3 (full output)");

  // scaling option
  int scale = 16;
  OCTAVE_HIGHS_GET_INT_PARAM ("scale", scale);
  if (scale < 0 || scale > 128)
    error ("__highs__: PARAM.scale must either be 128 (automatic selection of scaling options), or a bitwise or of: 1 (geometric mean scaling), 16 (equilibration scaling), 32 (round scale factors to power of two), 64 (skip if problem is well scaled");

  // Dual simplex option
  par.dual = 1;
  OCTAVE_HIGHS_GET_INT_PARAM ("dual", par.dual);
  if (par.dual < 1 || par.dual > 3)
    error ("__highs__: PARAM.dual must be 1 (use two-phase primal simplex [default]) or 2 (use two-phase dual simplex) or 3 (use two-phase dual simplex, and if it fails, switch to the primal simplex)");

  // Pricing option
  par.price = 34;
  OCTAVE_HIGHS_GET_INT_PARAM ("price", par.price);
  if (par.price != 17 && par.price != 34)
    error ("__highs__: PARAM.price must be 17 (textbook pricing) or 34 (steepest edge pricing [default])");

  // Simplex iterations limit
  par.itlim = std::numeric_limits<int>::max ();
  OCTAVE_HIGHS_GET_INT_PARAM ("itlim", par.itlim);

  // Output frequency, in iterations
  par.outfrq = 200;
  OCTAVE_HIGHS_GET_INT_PARAM ("outfrq", par.outfrq);

  // Branching heuristic option
  par.branch = 4;
  OCTAVE_HIGHS_GET_INT_PARAM ("branch", par.branch);
  if (par.branch < 1 || par.branch > 5)
    error ("__highs__: PARAM.branch must be 1 (first fractional variable) or 2 (last fractional variable) or 3 (most fractional variable) or 4 (heuristic by Driebeck and Tomlin [default]) or 5 (hybrid pseudocost heuristic)");

  // Backtracking heuristic option
  par.btrack = 4;
  OCTAVE_HIGHS_GET_INT_PARAM ("btrack", par.btrack);
  if (par.btrack < 1 || par.btrack > 4)
    error ("__highs__: PARAM.btrack must be 1 (depth first search) or 2 (breadth first search) or 3 (best local bound) or 4 (best projection heuristic [default]");

  // Presolver option
  par.presol = 1;
  OCTAVE_HIGHS_GET_INT_PARAM ("presol", par.presol);
  if (par.presol < 0 || par.presol > 1)
    error ("__highs__: PARAM.presol must be 0 (do NOT use LP presolver) or 1 (use LP presolver [default])");

  // LPsolver option
  int lpsolver = 1;
  OCTAVE_HIGHS_GET_INT_PARAM ("lpsolver", lpsolver);
  if (lpsolver < 1 || lpsolver > 2)
    error ("__highs__: PARAM.lpsolver must be 1 (simplex method) or 2 (interior point method)");

  // Ratio test option
  par.rtest = 34;
  OCTAVE_HIGHS_GET_INT_PARAM ("rtest", par.rtest);
  if (par.rtest != 17 && par.rtest != 34)
    error ("__highs__: PARAM.rtest must be 17 (standard ratio test) or 34 (Harris' two-pass ratio test [default])");

  par.tmlim = std::numeric_limits<int>::max ();
  OCTAVE_HIGHS_GET_INT_PARAM ("tmlim", par.tmlim);

  par.outdly = 0;
  OCTAVE_HIGHS_GET_INT_PARAM ("outdly", par.outdly);

  // Save option
  int save_pb = 0;
  OCTAVE_HIGHS_GET_INT_PARAM ("save", save_pb);
  save_pb = save_pb != 0;

  // Real parameters

  // Relative tolerance used to check if the current basic solution
  // is primal feasible
  par.tolbnd = 1e-7;
  OCTAVE_HIGHS_GET_REAL_PARAM ("tolbnd", par.tolbnd);

  // Absolute tolerance used to check if the current basic solution
  // is dual feasible
  par.toldj = 1e-7;
  OCTAVE_HIGHS_GET_REAL_PARAM ("toldj", par.toldj);

  // Relative tolerance used to choose eligible pivotal elements of
  //  the simplex table in the ratio test
  par.tolpiv = 1e-10;
  OCTAVE_HIGHS_GET_REAL_PARAM ("tolpiv", par.tolpiv);

  par.objll = -std::numeric_limits<double>::max ();
  OCTAVE_HIGHS_GET_REAL_PARAM ("objll", par.objll);

  par.objul = std::numeric_limits<double>::max ();
  OCTAVE_HIGHS_GET_REAL_PARAM ("objul", par.objul);

  par.tolint = 1e-5;
  OCTAVE_HIGHS_GET_REAL_PARAM ("tolint", par.tolint);

  par.tolobj = 1e-7;
  OCTAVE_HIGHS_GET_REAL_PARAM ("tolobj", par.tolobj);

  // Assign pointers to the output parameters
  ColumnVector xmin (mrowsc, octave_NA);
  double fmin = octave_NA;
  ColumnVector lambda (mrowsA, octave_NA);
  ColumnVector redcosts (mrowsc, octave_NA);

  double time = 0.0;
  int status = -1;

  int errnum = highs (sense, mrowsc, mrowsA, c, nz, rn.fortran_vec (),
                     cn.fortran_vec (), a.fortran_vec (), b, ctype,
                     freeLB.fortran_vec (), lb, freeUB.fortran_vec (),
                     ub, vartype.fortran_vec (), isMIP, lpsolver,
                     save_pb, scale, par, xmin.fortran_vec (), fmin,
                     status, lambda.fortran_vec (),
                     redcosts.fortran_vec (), time);

  octave_scalar_map extra;

  if (! isMIP)
    {
      extra.assign ("lambda", lambda);
      extra.assign ("redcosts", redcosts);
    }

  extra.assign ("time", time);
  extra.assign ("status", status);

  return ovl (xmin, fmin, errnum, extra);
}

/*
## No test needed for internal helper function.
%!assert (1)
*/

OCTAVE_NAMESPACE_END
