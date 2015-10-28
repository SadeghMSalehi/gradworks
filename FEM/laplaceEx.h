//
//  laplaceEx.h
//  kfem
//
//  Created by Joowhi Lee on 9/17/15.
//
//

#ifndef __kfem__laplaceEx__
#define __kfem__laplaceEx__

#include <stdio.h>
#include "kfem.h"

class LaplaceEx {
public:
	LaplaceEx(): fe(1), dof_handler(triangulation) {};
	
	void run();
	void make_grid();
	void setup_system();
	void assemble_system();
	void solve();
	void output_results() const;

	dealii::Triangulation<2> triangulation;
	dealii::FE_Q<2> fe;
	dealii::DoFHandler<2> dof_handler;
	
	dealii::SparsityPattern sparsity_pattern;
	dealii::SparseMatrix<double> system_matrix;
	dealii::Vector<double> solution;
	dealii::Vector<double> system_rhs;
};

#endif /* defined(__kfem__laplaceEx__) */
