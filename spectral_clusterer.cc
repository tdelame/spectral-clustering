/*  Created on: Jan 8, 2016
 *      Author: T. Delame (tdelame@gmail.com)
 */
# include <spectral_clusterer.h>
# include <omp_kmean.h>
# include <color_lut.h>

# include <spectra/SymEigsSolver.h>
# include <spectra/MatOp/SparseGenMatProd.h>

# include <algorithm>
# include <list>
# include <fstream>
# include <iomanip>

# include <nanoflann.h>

# include <omp.h>
BEGIN_PROJECT_NAMESPACE

  /*******************
   * Nanoflann setup *
   *******************/
  template< class BBOX >
  static void merge_bbox( BBOX& bb, real* point, size_t dimension )
  {
    for( size_t d = 0; d < dimension; ++ d )
      {
        bb[ d ].low = std::min( bb[ d ].low, point[ d ] );
        bb[ d ].high = std::max( bb[ d ].high, point[ d ] );
      }
  }

  template< class BBOX >
  static void merge_bbox( BBOX& out, BBOX& in, size_t dimension )
  {
    for( size_t d = 0; d < dimension; ++ d )
      {
        out[ d ].low = std::min( out[ d ].low, in[ d ].low );
        out[ d ].high = std::max( out[ d ].high, in[ d ].high );
      }
  }

  struct nanoflann_dataset
  {
    nanoflann_dataset( real* input, size_t n, size_t k ) :
      pts{ input }, size{ n }, dimension{k}
    {}

    real* pts;
    const size_t size;
    const size_t dimension;


    // Must return the number of data points
    inline size_t kdtree_get_point_count() const
    {
      return size;
    }

    inline real kdtree_distance(const real* a, const size_t b_idx, size_t size) const
    {
      real result = 0;
      real* b = pts + b_idx * size;
      for( size_t i = 0; i < size; ++ i, ++ b )
        {
          real temp = a[i] - *b;
          temp *= temp;
          result += temp;
        }
      return result;
    }


    inline real kdtree_get_pt(const size_t idx, int component) const
    {
      return pts[ idx * dimension + component ];
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& bb) const
    {
      for( size_t d = 0; d < dimension; ++ d )
        {
          bb[ d ].low = bb[ d ].high = pts[ d ];
        }

      # pragma omp parallel
      {
        BBOX thread_bb = bb;
        # pragma omp for
        for( size_t i = 1; i < size; ++ i )
          {
            merge_bbox( thread_bb, pts + i * dimension, dimension );
          }
        # pragma omp critical
        merge_bbox( bb, thread_bb, dimension );
      }
      return true;
    }
  };

  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor< real, nanoflann_dataset >,
    nanoflann_dataset
    > index_type;

  /******************************
   * Eigen triplet manipulation *
   ******************************/
  typedef Eigen::Triplet<real> T;

  bool triplet_order( const T& a, const T& b )
  {
    if( a.col() == b.col() )
      return a.row() < b.row();
    return a.col() < b.col();
  }

  bool triplet_unique( const T& a, const T& b )
  {
    return a.col() == b.col() && a.row() == b.row();
  }

  /***
   * Debug macros
   */

  static void print_2d_graph( Eigen::SparseMatrix<real>& weight, real* input_points )
  {
    auto source = input_points;

    auto minx = source[0];
    auto maxx = minx;
    auto miny = source[1];
    auto maxy = miny;
    source += 2;

    for( size_t i = 1; i < weight.cols(); ++ i, source += 2 )
      {
        auto x = source[0];
        if( minx > x )
          minx = x;
        else if( maxx < x )
          maxx = x;

        auto y = source[1];
        if( miny > y )
          miny = y;
        else if ( maxy < y )
          maxy = y;
      }
    auto scale = 800.0 / std::max( maxx - minx, maxy - miny );


    std::ofstream output( "output/similarity_graph.ps");
    output
      << "%!PS-Adobe-3.0 EPSF-3.0\n"
      << "%%BoundingBox: -5 -5 800 800\n"
      << "/l {rlineto} def\n"
      << "/m {rmoveto} def\n"
      << "/c { .25 sub exch .25 sub exch 1.5 0 360 arc fill } def\n"
      << "/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
      << " gsave 1 setgray fill grestore gsave 3 setlinewidth"
      << " 1 setgray stroke grestore 0 setgray stroke } def\n";

    output << "0.2 0.5 0.8 setrgbcolor\n";
    for( size_t i = 0; i < weight.cols(); ++ i )
      {
        output << input_points[ i << 1 ] * scale << " " << input_points[ i * 2 + 1 ] * scale << " c\n";
      }

    output << "0.2 0.8 0.1 setrgbcolor\n";
    for( int k = 0; k < weight.outerSize(); ++ k )
       for( Eigen::SparseMatrix<real>::InnerIterator it(weight,k); it; ++ it )
         {
           if( it.row() > it.col() ) continue;
           output << "newpath " << input_points[ it.row() * 2 ] * scale << " " << input_points[ it.row() * 2 + 1] * scale << " moveto "
                  << input_points[ it.col() * 2 ] * scale << " " << input_points[ it.col() * 2 + 1] * scale << " lineto stroke\n";
         }
    output << "showpage\n";
    output.close();
  }

  static void print_kmean_on_2d_eigenvectors( std::vector<real>& eigenvectors,  omp_kmean& alg )
  {
    auto source = eigenvectors.data();

    auto minx = source[0];
    auto maxx = minx;
    auto miny = source[1];
    auto maxy = miny;

    for( size_t i = 1; i < alg.get_number_of_points(); ++ i, source += 2 )
      {
        auto x = source[0];
        if( minx > x )
          minx = x;
        else if( maxx < x )
          maxx = x;

        auto y = source[1];
        if( miny > y )
          miny = y;
        else if ( maxy < y )
          maxy = y;
      }
    auto scale = 800.0 / std::max( maxx - minx, maxy - miny );
    std::ofstream output( "output/eigenvector_kmean.ps");
    output
      << "%!PS-Adobe-3.0 EPSF-3.0\n"
      << "%%BoundingBox: " << minx * scale<< " " << miny* scale << " " << maxx* scale << " " << maxy * scale<< "\n"
      << "/l {rlineto} def\n"
      << "/m {rmoveto} def\n"
      << "/c { .25 sub exch .25 sub exch 1.5 0 360 arc fill } def\n"
      << "/centroid { .25 sub exch .25 sub exch 2.5 0 360 arc fill } def\n"
      << "/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
      << " gsave 1 setgray fill grestore gsave 3 setlinewidth"
      << " 1 setgray stroke grestore 0 setgray stroke } def\n";


    for( size_t i = 0; i < alg.get_number_of_points(); ++ i )
      {
        if( alg.get_group( i ) )
          output << "0.2 0.4 0.8 setrgbcolor\n";
        else output << "0.2 0.8 0.4 setrgbcolor\n";
        output << eigenvectors[ i * 2 ]* scale << " " << eigenvectors[ i * 2 + 1]* scale << " c\n";
      }
    for( size_t k = 0; k < alg.get_number_of_clusters(); ++ k )
      {
        output << "0.8 0.4 0.2 setrgbcolor\n";
        output << alg.get_centroid( k )[0]* scale << " " << alg.get_centroid( k )[1]* scale << " centroid\n";
        std::cout << "centroid " << k << ": " << alg.get_centroid( k )[0] << " " << alg.get_centroid( k )[1] << std::endl;
      }

    output << "showpage\n";
    output.close();
  }

  static void print_weight_matrix( Eigen::SparseMatrix<real>& weights )
  {
    typedef unsigned char uchar;
    auto max_weight = real(0);
    auto min_weight = REAL_MAX;
    for( int k = 0; k < weights.outerSize(); ++ k )
      for( Eigen::SparseMatrix<real>::InnerIterator it(weights,k); it; ++ it )
        {
          max_weight = std::max( max_weight, it.value() );
          min_weight = std::min( min_weight, it.value() );
        }

    uchar* pixels = new uchar[ weights.cols() * weights.cols() * 3 ];
    std::memset( pixels, uchar(255), weights.cols() * weights.cols() * 3);
    real color[3];
    for( int k = 0; k < weights.outerSize(); ++ k )
      for( Eigen::SparseMatrix<real>::InnerIterator it(weights,k); it; ++ it )
        {
          set_color( it.value(), min_weight, max_weight, color );
          auto dest = pixels + (it.row() * weights.cols() + it.col()) * 3;
          dest[0] = uchar( color[0] * 255.0 );
          dest[1] = uchar( color[1] * 255.0 );
          dest[2] = uchar( color[2] * 255.0 );
        }
    std::ofstream ppm ("output/weight_matrix.ppm", std::ios::binary);
    ppm << "P6\n" << int(weights.cols()) << " " << int(weights.cols()) << "\n255\n";
    auto dest = pixels;
    for( size_t i = weights.cols() * weights.cols() * 3; i; --i, ++ dest )
      {
        ppm << *dest;
      }
    delete[] pixels;
    ppm.close();
  }

  spectral_clusterer::parameters::parameters()
  : m_graph_mode{ DISTANCE_GRAPH }, m_laplacian_mode{ LAPLACIAN },
    m_distance_threshold{ 1e-2 }, m_similarity_threshold{ 1e-6 },
    m_number_of_neighbors{ 5 }, m_number_of_clusters{ 3 }
  {}

  spectral_clusterer::spectral_clusterer(
      real* input, size_t number_of_points,
      size_t dimension_of_points,
      const parameters& params )
  : m_input{ input }, m_n{ number_of_points }, m_d{ dimension_of_points },
    m_params{ params }, m_weights( m_n, m_n ), diagonal( m_n, 0 ), m_diagonal( m_n, m_n ),
    m_eigenvectors( params.m_number_of_clusters * m_n, 0 ), m_duration{ omp_get_wtime() }
  {
    build_weight_matrix();
    extract_eigen_vectors();
    kmean_in_coordinate_space();
    m_duration = omp_get_wtime() - m_duration;
  }

  real
  spectral_clusterer::get_execution_duration() const
  {
    return m_duration;
  }

  size_t
  spectral_clusterer::get_group( size_t index ) const
  {
    return m_groups[ index ];
  }

  size_t
  spectral_clusterer::get_number_of_clusters() const
  {
    return m_params.m_number_of_clusters;
  }

  void
  spectral_clusterer::build_weight_matrix()
  {
    std::list< T > triplets;

    nanoflann_dataset dataset( m_input, m_n, m_d );
    index_type index( m_d, dataset, nanoflann::KDTreeSingleIndexAdaptorParams{32} );
    index.buildIndex();

    switch( m_params.m_graph_mode )
    {
      case parameters::DISTANCE_GRAPH:
        {
          nanoflann::SearchParams index_search( -1, 0.0, false );
          real squared_radius = m_params.m_distance_threshold * m_params.m_distance_threshold;
          real inv_two_sigma_squared = 1.0 / ( squared_radius );
          # pragma omp parallel
          {
            std::list< T > thread_triplets;
            std::vector< std::pair< size_t, real > > results;
            results.resize( m_n );
            # pragma omp for
            for( size_t i = 0; i < m_n; ++ i )
              {
                size_t count = index.radiusSearch( m_input + i * m_d, squared_radius, results, index_search );
                for( size_t k = 0; k < count; ++ k )
                  {
                    auto& pair = results[k];
                    if( pair.first < i )
                      {
                        auto similarity = std::exp( -pair.second * inv_two_sigma_squared );
                        thread_triplets.push_back( T( i, pair.first, similarity ) );
                        thread_triplets.push_back( T( pair.first, i, similarity ) );
                      }
                  }
              }
            # pragma omp critical
            triplets.splice( triplets.end(), thread_triplets );
          }
          m_weights.setFromTriplets( triplets.begin(), triplets.end() );
          break;
        }
      case parameters::KNN_GRAPH:
        {
          const size_t k = m_params.m_number_of_neighbors;
          # pragma omp parallel
          {
            std::list< T > thread_triplets;
            std::vector< size_t > indices( k, 0 );
            std::vector< real > sdistances( k, 0 );

            # pragma omp for schedule(dynamic)
            for( size_t i = 0; i < m_n; ++ i )
              {
                index.knnSearch(
                    m_input + i * m_d,
                    m_params.m_number_of_neighbors,
                    indices.data(),
                    sdistances.data() );

                for( size_t j = 0; j < k; ++ j )
                  {
                    if( indices[j] == i ) continue;
                    auto similarity = 1.0;
                    if( i < indices[j] )
                      thread_triplets.push_back( T( i, indices[j], similarity ) );
                    else
                      thread_triplets.push_back( T( indices[j], i, similarity ) );
                  }
              }

            thread_triplets.sort( triplet_order );
            # pragma omp critical
            triplets.merge( thread_triplets, triplet_order );
          }
          triplets.unique( triplet_unique );
          std::list< T > symmetry;
          for( auto& t : triplets )
            symmetry.push_back( T( t.col(), t.row(), t.value() ) );
          triplets.splice( triplets.end(), symmetry );
          m_weights.setFromTriplets( triplets.begin(), triplets.end() );
          break;
        }
      case parameters::SYMMETRIC_KNN_GRAPH:
        {
          const size_t k = m_params.m_number_of_neighbors;
          # pragma omp parallel
          {
            std::list< T > thread_triplets;
            std::vector< size_t > indices( k, 0 );
            std::vector< real > sdistances( k, 0 );

            # pragma omp for schedule(dynamic)
            for( size_t i = 0; i < m_n; ++ i )
              {
                index.knnSearch(
                    m_input + i * m_d,
                    m_params.m_number_of_neighbors,
                    indices.data(),
                    sdistances.data() );

                for( size_t j = 0; j < k; ++ j )
                  {
                    if( indices[j] == i ) continue;
                    auto similarity = 1.0;
                    if( i < indices[j] )
                      thread_triplets.push_back( T( i, indices[j], similarity ) );
                    else
                      thread_triplets.push_back( T( indices[j], i, similarity ) );
                  }
              }

            thread_triplets.sort( triplet_order );
            # pragma omp critical
            triplets.merge( thread_triplets, triplet_order );
          }
          auto it = triplets.begin(), end = triplets.end();
          while( it != end )
            {
              auto next = it;
              ++next;

              // check if next is the duplicate of it
              if( next != end )
                {
                  // they are duplicates, thus this connection is symmetrical
                  if( triplet_unique( *it, *next ) )
                    {
                      // instead of removing next and then adding the symmetrical value
                      // just update next to represent the symmetrical value
                      *next = T( it->col(), it->row(), it->value() );
                      it = next;
                      ++it;
                    }
                  else
                    {
                      it = triplets.erase( it );
                    }
                }
              // it was the last one, thus has no duplicate
              else
                it = triplets.erase( it );
            }
          m_weights.setFromTriplets( triplets.begin(), triplets.end() );
          break;
        }
      }

    for( size_t i = 0; i < m_n; ++ i )
      {
        diagonal[ i ] = m_weights.col(i).sum();
      }
# ifdef DEBUG
    print_weight_matrix(m_weights);
    if( m_d == 2 )
      {
        print_2d_graph( m_weights, m_input );
      }
# endif
  }

  void
  spectral_clusterer::compute_first_eigenvectors( Eigen::SparseMatrix<real>& laplacian )
  {
    /* A parameter that controls the convergence of the spectra algorithm.
     * A higher value means a better chance to obtain the results before the
     * maximum iteration number is reached, but also a higher memory consumption.
     * According to the implementation, this number should be at least equal to
     * the double of the number of eigen values requested for the algorithm to be
     * efficient.*/
    int ncv = 4 * m_params.m_number_of_clusters;
    if( ncv > m_n )
      ncv = std::min( int(2 * m_params.m_number_of_clusters), int(m_n));

    Spectra::SparseGenMatProd<real> op(laplacian);
    Spectra::SymEigsSolver< real, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<real> > eigs(
        &op, m_params.m_number_of_clusters, ncv);
    eigs.init();

    int iterations = 1000;
    while( true )
      {
        int nconv = eigs.compute( iterations );
        if( nconv == m_params.m_number_of_clusters )
          {
            break;
          }
        iterations *= 10;
        LOG( info, "Only " << nconv <<" eigen value(s) found. Raising the max number of iterations to " << iterations );
      }

    auto temp = eigs.eigenvectors();
    for( size_t i = 0; i < m_n; ++ i )
      {
        for( size_t k = 0; k < m_params.m_number_of_clusters; ++ k )
          {
            m_eigenvectors[ i * m_params.m_number_of_clusters + k ] = temp( i, k );
          }
      }
  }

  void
  spectral_clusterer::extract_eigen_vectors()
  {
    if( m_params.m_laplacian_mode == parameters::LAPLACIAN )
      {
        Eigen::SparseMatrix<real> laplacian = - m_weights;
        for( size_t i = 0; i < m_n; ++ i )
          laplacian.coeffRef( i, i ) = diagonal[ i ];
        compute_first_eigenvectors( laplacian );
      }
    else
      {
        Eigen::SparseMatrix<real> laplacian( m_n, m_n );
        std::list< T > triplets;
        std::vector< real > inv_sqrt_diagonal = diagonal;
        for( size_t i = 0; i < m_n; ++ i )
          {
            triplets.push_back( T( i, i, real(1.0 ) ) );
            inv_sqrt_diagonal[i] = real(1.0)/ std::sqrt( inv_sqrt_diagonal[i] );
          }

        for( int k = 0; k < m_weights.outerSize(); ++ k )
          {
            for( Eigen::SparseMatrix<real>::InnerIterator it(m_weights,k); it; ++ it )
              {
                triplets.push_back(
                    T(
                      it.row(),
                      it.col(),
                      -inv_sqrt_diagonal[it.row()] * it.value() * inv_sqrt_diagonal[it.col()] ) );
              }
          }
        laplacian.setFromTriplets( triplets.begin(), triplets.end() );
        compute_first_eigenvectors( laplacian );

        if( m_params.m_laplacian_mode == parameters::RANDOM_WALK_LAPLACIAN )
          {
            // instead of using yet another library to compute the generalized eigen problem:
            // (D - W) * u = lambda * diagonal * u (1),
            // I use the relationship between the solutions to (1) and the solutions to:
            // (I - D^(-1/2) * W * D^(-1/2)) * u = lambda * u
            auto dst = m_eigenvectors.data();
            for( size_t i = 0; i < m_n; ++ i )
              {
                for( size_t k = 0; k < m_params.m_number_of_clusters; ++k, ++dst )
                  {
                    *dst *= inv_sqrt_diagonal[ i ];
                  }
              }
          }
        else
          {
            // this method requires an additional step: normalization of rows
            for( size_t i = 0; i < m_n; ++ i )
              {
                auto factor = real(0.0);
                for( size_t j = 0; j < m_params.m_number_of_clusters; ++ j )
                  {
                    auto val = m_eigenvectors[ i * m_params.m_number_of_clusters + j ];
                    val *= val;
                    factor += val;
                  }
                factor = real(1.0) / std::sqrt( factor );
                for( size_t j = 0; j < m_params.m_number_of_clusters; ++ j )
                  {
                    m_eigenvectors[ i * m_params.m_number_of_clusters + j ] *= factor;
                  }
              }
          }
      }
  }

  void
  spectral_clusterer::kmean_in_coordinate_space()
  {
    omp_kmean alg(
        m_eigenvectors.data(),
        m_n,
        m_params.m_number_of_clusters,
        m_params.m_number_of_clusters,
        m_groups );
    alg.execute();
# ifdef DEBUG
    print_kmean_on_2d_eigenvectors( m_eigenvectors, alg );
# endif

  }

END_PROJECT_NAMESPACE
