/*  Created on: Dec 28, 2015
 *      Author: T. Delame (tdelame@gmail.com)
 */
# include <project.h>
# include <color_lut.h>
# include <spectral_clusterer.h>

# include <valgrind/callgrind.h>

# include <fstream>
# include <iomanip>

BEGIN_PROJECT_NAMESPACE

  class data_interface {
  public:
    data_interface( size_t n, real diagonal )
      : m_size{ n }, m_diagonal{ diagonal }
    {}

    virtual ~data_interface()
    {}

    operator real*()
    {
      return get_data_pointer();
    }
    size_t get_size() const
    {
      return m_size;
    }
    real get_diagonal() const
    {
      return m_diagonal;
    }

  protected:
    virtual real* get_data_pointer() = 0;
    const size_t m_size;
    const real m_diagonal;
  };

  static
  void fill_circle(
      real* destination,
      size_t number_of_points,
      real center[2],
      real radii[2] )
  {
    static real two_pi = 2.0 * M_PI;
    real thickness = radii[1] - radii[0];
    for( size_t i = 0; i < number_of_points; ++ i, destination += 2 )
      {
        auto angle = unit_random() * two_pi;
        auto r = radii[0] + unit_random() * thickness;
        destination[ 0 ] = center[0] + r * std::cos(angle);
        destination[ 1 ] = center[1] + r * std::sin(angle);
      }
  }

  class disks_data
    : public data_interface {
  public:
    disks_data( size_t n, real diagonal )
      : data_interface{ n, diagonal },
        m_data{ new real[ n * 2 ] }
    {
      auto n1 = n >> 1;
      real center[2] = { m_diagonal / 2, m_diagonal / 2 };
      real radii[2] = { m_diagonal * 0.2, m_diagonal * 0.3 };
      fill_circle( m_data, n1, center, radii );

      auto n2 = n - n1;
      radii[0] = 0.4 * m_diagonal;
      radii[1] = 0.45 * m_diagonal;
      fill_circle( m_data + 2 * n1, n2, center, radii );
    }

    ~disks_data()
    {
      delete[] m_data;
    }

  protected:
    real* get_data_pointer() override
    {
      return m_data;
    }

  private:
    real* m_data;
  };

  static void
  save_results(
      data_interface& data,
      spectral_clusterer& alg,
      const std::string& name )
  {
    const auto nb_clusters = alg.get_number_of_clusters();
    std::vector< real > colors( 3 * data.get_size(), 0 );
    for( size_t i = 0; i < nb_clusters; ++ i )
      {
        set_color( real(i), real(0), real(nb_clusters), colors.data() + 3 * i );
      }

    const auto scale = real(800) / ( data.get_diagonal() );

    std::ofstream output( name );

    output
      << "%!PS-Adobe-3.0 EPSF-3.0\n"
      << "%%BoundingBox: -5 -5 800 800\n"
      << "/l {rlineto} def\n"
      << "/m {rmoveto} def\n"
      << "/c { .25 sub exch .25 sub exch 2.5 0 360 arc fill } def\n"
      << "/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
      << " gsave 1 setgray fill grestore gsave 3 setlinewidth"
      << " 1 setgray stroke grestore 0 setgray stroke } def\n";

    for( size_t cluster_id = 0; cluster_id < nb_clusters; ++ cluster_id )
      {
        output << colors[ 3 * cluster_id ] << " "
            << colors[ 3 * cluster_id + 1 ] << " "
            << colors[ 3 * cluster_id + 2 ] << " setrgbcolor\n";
        real* source = data;
        for( size_t i = 0; i < data.get_size(); ++ i, source += 2 )
          {
            if( alg.get_group( i ) == cluster_id )
              {
                output << scale * source[0] << " " << scale * source[ 1 ] << " c\n";
              }
          }
      }
     output << "showpage\n%%EOF";
     output.close();
  }

  static int
  do_execute( int argc, char* argv[] )
  {
    (void)argc;
    (void)argv;

    disks_data data( 4000, 40 );

    spectral_clusterer::parameters params;
    params.m_graph_mode = spectral_clusterer::parameters::DISTANCE_GRAPH;
    params.m_number_of_clusters = 2;
    params.m_distance_threshold = data.get_diagonal() * 0.1;
    params.m_number_of_neighbors = 8;

//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/1_disk_DG_L_optimal.ps");
//      LOG( info, "DG_L_optimal took " << alg.get_execution_duration() );
//    }
//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::SYMMETRIC_LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/2_disk_DG_LS_optimal.ps");
//      LOG( info, "DG_LS_optimal took " << alg.get_execution_duration() );
//    }
//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::RANDOM_WALK_LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/3_disk_DG_LRW_optimal.ps");
//      LOG( info, "DG_LRW_optimal took " << alg.get_execution_duration() );
//    }
//
//    params.m_graph_mode = spectral_clusterer::parameters::KNN_GRAPH;
//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/4_disk_KG_L_5_param.ps");
//      LOG( info, "KG_L_5 took " << alg.get_execution_duration() );
//    }
//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::SYMMETRIC_LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/5_disk_KG_LS_5_param.ps");
//      LOG( info, "KG_LS_5 took " << alg.get_execution_duration() );
//    }
//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::RANDOM_WALK_LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/6_disk_KG_LRW_5_param.ps");
//      LOG( info, "KG_LRW_5 took " << alg.get_execution_duration() );
//    }
//
//    params.m_graph_mode = spectral_clusterer::parameters::SYMMETRIC_KNN_GRAPH;
//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/7_disk_SKG_L_5_param.ps");
//      LOG( info, "SKG_L_5 took " << alg.get_execution_duration() );
//    }
//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::SYMMETRIC_LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/8_disk_SKG_LS_5_param.ps");
//      LOG( info, "SKG_LS_5 took " << alg.get_execution_duration() );
//    }
//    {
//      params.m_laplacian_mode = spectral_clusterer::parameters::RANDOM_WALK_LAPLACIAN;
//      spectral_clusterer alg( data, data.get_size(), 2, params );
//      save_results( data, alg, "output/9_disk_SKG_LRW_5_param.ps");
//      LOG( info, "SKG_LRW_5 took " << alg.get_execution_duration() );
//    }

    {
      params.m_laplacian_mode = spectral_clusterer::parameters::RANDOM_WALK_LAPLACIAN;

      const auto internal = real(1.0);
      const auto external = real(0.8);
      const auto bridge   = real(0.2);

      weight_matrix weights;
      weights.m_size = 16;


      weights.m_coeffs = {
          {0,1,internal}, {1,2,internal}, {2,3,internal}, {3,4,internal},
          {4,5,internal}, {5,6,internal}, {6,7,internal}, {7,8,internal},
          {1,0,internal}, {2,1,internal}, {3,2,internal}, {4,3,internal},
          {5,4,internal}, {6,5,internal}, {7,6,internal}, {8,7,internal},

          {0,9,bridge}, {9,0,bridge},

          {8,9,external}, {9,10,external}, {10,11,external}, {11,12,external},
          {12,13,external}, {13,14,external}, {14,15,external}, {15,8,external},
          {9,8,external}, {10,9,external}, {11,10,external}, {12,11,external},
          {13,12,external}, {14,13,external}, {15,14,external}, {8,15,external}
      };

      spectral_clusterer alg( weights, params );
      LOG( info, "Graph LRW took " << alg.get_execution_duration() );
      for( size_t i = 0; i < weights.m_size; ++ i )
        {
          std::cout << "node " << std::setw(3) << i << ": " << alg.get_group( i ) << "\n";
        }
    }

    return EXIT_SUCCESS;
  }

END_PROJECT_NAMESPACE

int main( int argc, char* argv[] )
{
  return PROJECT_NAMESPACE::do_execute( argc, argv );
}
