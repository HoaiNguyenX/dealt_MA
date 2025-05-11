


// strong residual error estimator
template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    minimal_residual_error_estimate(
      const std::vector< unsigned int >&  degrees,
      const Vector<double>&               solution,
      std::map< cell_iterator, double >&  residuals
    ) const 
  {
    Assert(this->is_bezier_mesh, ExcInvalidState());

    TSValues<dim, spacedim> ts_values(
      this, degrees,
      update_values |
      update_gradients |
      update_hessians |
      update_quadrature_points |
      update_JxW_values
    );
    TSFaceValues<dim, spacedim> face_values(
      this, degrees,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    const unsigned int nvc = GeometryInfo<2>::vertices_per_cell;
    const unsigned int nvf = GeometryInfo<dim>::vertices_per_face;

    std::map< unsigned int, double > face_residuals;
    for (const auto& cell : this -> active_cell_iterators()) {
      std::vector< unsigned int > local_dof_indices = get_IEN_array(cell);   

      // Compute the jumps of faces at C0 continuity
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        // Get orientation of current face:
        const auto& face = cell -> face(f);
        const Point<dimension>& c = (-1.) * 
                                    face -> vertex(0) +
                                    face -> vertex(nvf - 1);
        unsigned int orientation = 0;
        for ( ; orientation < dimension && c(orientation) != 0; orientation++);

        if (face -> has_children() ) {
          // This is merely a special case, as refined faces are considered
          // in-active and hence their index is not represented in mof
          if (mof.at(face -> child(0) -> index()) == p[orientation]) {
            // C0 Edge detected:
            std::cout << "C0 edge detected 1" << std::endl;

            face_values.reinit(cell, f);
            const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
            for (const unsigned int q_index : face_values.quadrature_point_indices()){
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                        face_values.shape_grad(i, q_index) *
                        face_values.normal_vector(q_index);

              const unsigned int ch = q_index > ppf ? 1 : 0;
              face_residuals[face -> child(ch) -> index()] +=
                      0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // otherwise it is atleast a C1 edge, and cannot be at the boundary
        } else if (mof.at(face -> index()) == p[orientation]){
          // C0 edge detected
          face_values.reinit(cell, f);

          // Since this face is not refined, we can simply compute
          // the jump terms along it and store the values at the
          // corresponding place in face_residuals
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = 0;
            for (const unsigned int i : face_values.dof_indices())
              g += solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += 0.25 * g * g * face_values.JxW(q_index);
          } // for ( q_index )
          // We dont have neumann bc here
        } // if ( ... )
      } // for ( face )
    } // for ( cell )
    
    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals[cell];
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);

      // get old solution gradient
      std::vector< Tensor<1, spacedim> > grad_u(ts_values.n_quadrature_points_per_cell());
      ts_values.get_function_gradients(solution, grad_u);

      // get old solution hessian
       std::vector< Tensor<2, spacedim> > hessian(ts_values.n_quadrature_points_per_cell());
       ts_values.get_function_hessians(solution, hessian);
      
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){

        double first_term = 0., 
               second_term = 0.;
        for (const unsigned int i : ts_values.dof_indices() )
        {
          first_term += (1 + grad_u[q_index].norm_square())
                         * laplace_u;
          second_term += grad_u[q_index] * hessian_u[q_index]
                        * grad_u[q_index];
        } // for ( i )
        double tmp = first_term - second_term;
        local_residual += tmp * tmp * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                                cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;

      // And then also add the values from the faces:
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        if (cell -> face(f) -> has_children()){
          const double child0_width = 
             (cell -> face(f) -> child(0) -> vertex(0)).distance(
              cell -> face(f) -> child(0) -> vertex(nvf - 1));
          const double child1_width = 
             (cell -> face(f) -> child(1) -> vertex(0)).distance(
              cell -> face(f) -> child(1) -> vertex(nvf - 1));
          local_residual += child0_width *
                  face_residuals[cell -> face(f) -> child(0) -> index()];
          local_residual += child1_width *
                  face_residuals[cell -> face(f) -> child(1) -> index()];
        } else {
          const double face_width = 
              (cell -> face(f) -> vertex(0)).distance(
               cell -> face(f) -> vertex(nvf - 1));
          local_residual += face_width *
                  face_residuals[cell -> face(f) -> index()];
        }
      } // for ( f )

      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // minimal_residual_error_estimate() [1/2]









// weak residual error estimator
  template<int dim, int spacedim>
  void TS_TriangulationBase<dim, spacedim>::
    minimal_residual_error_estimate(
      const std::vector< unsigned int >&  degrees,
      const Vector<double>&               solution,
      std::map< cell_iterator, double >&  residuals
    ) const 
  {
    Assert(this->is_bezier_mesh, ExcInvalidState());

    TSValues<dim, spacedim> ts_values(
      this, degrees,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values
    );
    TSFaceValues<dim, spacedim> face_values(
      this, degrees,
      update_values |
      update_gradients |
      update_quadrature_points |
      update_JxW_values  |
      update_normal_vectors);
    const unsigned int nvc = GeometryInfo<2>::vertices_per_cell;
    const unsigned int nvf = GeometryInfo<dim>::vertices_per_face;

    std::map< unsigned int, double > face_residuals;
    for (const auto& cell : this -> active_cell_iterators()) {
      std::vector< unsigned int > local_dof_indices = get_IEN_array(cell);   

      // Compute the jumps of faces at C0 continuity
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        // Get orientation of current face:
        const auto& face = cell -> face(f);
        const Point<dimension>& c = (-1.) * 
                                    face -> vertex(0) +
                                    face -> vertex(nvf - 1);
        unsigned int orientation = 0;
        for ( ; orientation < dimension && c(orientation) != 0; orientation++);

        if (face -> has_children() ) {
          // This is merely a special case, as refined faces are considered
          // in-active and hence their index is not represented in mof
          if (mof.at(face -> child(0) -> index()) == p[orientation]) {
            // C0 Edge detected:
            std::cout << "C0 edge detected 1" << std::endl;

            face_values.reinit(cell, f);
            const unsigned int ppf = 2 * face_values.n_quadrature_points_per_face();
            for (const unsigned int q_index : face_values.quadrature_point_indices()){
              double g = 0;
              for (const unsigned int i : face_values.dof_indices())
                g += solution(local_dof_indices[i]) *
                        face_values.shape_grad(i, q_index) *
                        face_values.normal_vector(q_index);

              const unsigned int ch = q_index > ppf ? 1 : 0;
              face_residuals[face -> child(ch) -> index()] +=
                      0.25 * g * g * face_values.JxW(q_index);
            } // for ( q_index )
          } // otherwise it is atleast a C1 edge, and cannot be at the boundary
        } else if (mof.at(face -> index()) == p[orientation]){
          // C0 edge detected
          face_values.reinit(cell, f);

          // Since this face is not refined, we can simply compute
          // the jump terms along it and store the values at the
          // corresponding place in face_residuals
          for (const unsigned int q_index : face_values.quadrature_point_indices()){
            double g = 0;
            for (const unsigned int i : face_values.dof_indices())
              g += solution(local_dof_indices[i]) *
                      face_values.shape_grad(i, q_index) *
                      face_values.normal_vector(q_index);

            face_residuals[face -> index()] += 0.25 * g * g * face_values.JxW(q_index);
          } // for ( q_index )
          // We dont have neumann bc here
        } // if ( ... )
      } // for ( face )
    } // for ( cell )
    
    // After every Jump on every face is computed,
    // we compute the residuals on the cells
    for (const auto& cell : this -> active_cell_iterators()){
      double& local_residual = residuals[cell];
      const std::vector< unsigned int >& local_dof_indices = get_IEN_array(cell);
      ts_values.reinit(cell);

      // get old solution gradient
      std::vector< Tensor<1, spacedim> > sol_grad(ts_values.n_quadrature_points_per_cell());
      ts_values.get_function_gradients(solution, sol_grad);
      
      for ( const unsigned int q_index : ts_values.quadrature_point_indices()){

        // coefficient a_n depending on old solution gradient for readability
        const double coeff = 1. / std::sqrt(1 + sol_grad[q_index] * sol_grad[q_index]);
        double g = 0.;
        for (const unsigned int i : ts_values.dof_indices() ){
          double tmp = ts_values.shape_grad(i, q_index)
                            * coeff     
                            * sol_grad[q_index];
          g += tmp * tmp;            
        } // for ( i )

        local_residual += g * ts_values.JxW(q_index);
      } // for ( q_index )
      const double cell_width = (cell->vertex(0)).distance(
                                cell->vertex(nvc - 1));

      local_residual *= cell_width * cell_width;

      // And then also add the values from the faces:
      for (unsigned int f = 0;
              f < GeometryInfo<dimension>::faces_per_cell;
              f++) {
        if (cell -> face(f) -> has_children()){
          const double child0_width = 
             (cell -> face(f) -> child(0) -> vertex(0)).distance(
              cell -> face(f) -> child(0) -> vertex(nvf - 1));
          const double child1_width = 
             (cell -> face(f) -> child(1) -> vertex(0)).distance(
              cell -> face(f) -> child(1) -> vertex(nvf - 1));
          local_residual += child0_width *
                  face_residuals[cell -> face(f) -> child(0) -> index()];
          local_residual += child1_width *
                  face_residuals[cell -> face(f) -> child(1) -> index()];
        } else {
          const double face_width = 
              (cell -> face(f) -> vertex(0)).distance(
               cell -> face(f) -> vertex(nvf - 1));
          local_residual += face_width *
                  face_residuals[cell -> face(f) -> index()];
        }
      } // for ( f )

      local_residual = std::sqrt(local_residual);
    } // for ( cell )
  } // minimal_residual_error_estimate() [1/2]