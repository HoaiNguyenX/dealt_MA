

















/*
IPF_Data<3> Nonlinear_Benchmark::get_IPF_data() const
  {
    std::vector< std::vector< double > > kv;
    std::vector< Point<3> >              cps;
    std::vector< unsigned int >          deg;

    // define the knot vectors
    kv = std::vector< std::vector< double > >(3);
    kv[0] = {0, 0, 0, 1, 1, 1};
    kv[1] = {0, 0, 0, 1, 1, 1};
    kv[2] = {0, 0, 0, 1, 1, 1};

    // define the control points vector
    cps = std::vector< Point<3> >(8);
    cps[0] = Point<3>( 0.,  0.,  0.);
    cps[1] = Point<3>( 1.,  0.,  0.);
    cps[2] = Point<3>( 0.,  1.,  0.);
    cps[3] = Point<3>( 1.,  1.,  0.);
    cps[4] = Point<3>( 0.,  0.,  1.);
    cps[5] = Point<3>( 1.,  0.,  1.);
    cps[6] = Point<3>( 0.,  1.,  1.);
    cps[7] = Point<3>( 1.,  1.,  1.);

    // define the degrees
    deg = {1, 1, 1};

    std::vector< double >     weights(8, 1.);
    
    return IPF_Data<3>(cps, weights, kv, deg);;
  } // get_IPF_data 3-dim


  */