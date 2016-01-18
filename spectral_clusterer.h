/*  Created on: Jan 8, 2016
 *      Author: T. Delame (tdelame@gmail.com)
 */

# ifndef SPECTRAL_CLUSTERER_H_
# define SPECTRAL_CLUSTERER_H_
# include <project.h>
# include <vector>

# include <Eigen/Core>
# include <Eigen/Sparse>

BEGIN_PROJECT_NAMESPACE

  /**
   * Spectral clustering on n d-dimensional points. The algorithm manages 3
   * different similarity graph construction modes, all using the euclidean
   * distance between points:
   *
   *   - DISTANCE_GRAPH: The similarity between points is the function
   *   exp( - 0.5 * d^2 / sigma^2 ), where d is the distance between two
   *   points. The edges are weighted by this similarity.
   *
   *   - KNN_GRAPH: the edges of a point are formed with its k nearest neighbors and
   *   all have a weight of 1
   *
   *   - SYMMETRIC_KNN_GRAPH: same as KNN_GRAPH except an edge between points p1 and
   *   p2 is in the graph only if p1 is in the KNN of p2 and p2 is in the KNN of p1.
   *
   * The weight matrix of the graph W is composed of the weights w_ij which are non
   * negative for existing edges (and null for every couple of points (i,j) not
   * connected by an edge). The diagonal matrix D = (d_ii) is defined by
   * d_ii = sum_j w_ij.
   *
   * The algorithm also manages 3 kinds of laplacian matrix mode:
   *
   *   - LAPLACIAN: L = D - W and we solve L * u = lambda * u
   *
   *   - RANDOM_WALK_LAPLACIAN: L = D - W and we solve L * u = lambda * D * u
   *
   *   - SYMMETRIC_LAPLACIAN: L = I - D^(-1/2) * W * D^(-1/2) and we solve
   *   L * u = lambda * u
   *
   */
  class spectral_clusterer {
  public:
    struct parameters {
      enum graph_mode {
        DISTANCE_GRAPH,
        KNN_GRAPH,
        SYMMETRIC_KNN_GRAPH
      };

      enum laplacian_mode {
        LAPLACIAN,
        RANDOM_WALK_LAPLACIAN,
        SYMMETRIC_LAPLACIAN
      };

      graph_mode m_graph_mode;
      laplacian_mode m_laplacian_mode;
      real m_distance_threshold;
      real m_similarity_threshold;
      unsigned int m_number_of_neighbors;
      unsigned int m_number_of_clusters;

      parameters();
    };

    spectral_clusterer(
      real* input,
      size_t number_of_points,
      size_t dimension_of_points,
      const parameters& params = parameters{} );

    size_t
    get_group( size_t index ) const;

    size_t
    get_number_of_clusters() const;

    real
    get_execution_duration() const;

  private:
    void build_weight_matrix();
    void extract_eigen_vectors();
    void kmean_in_coordinate_space();

    void compute_first_eigenvectors( Eigen::SparseMatrix<real>& laplacian );

    real* m_input;
    const size_t m_n;
    const size_t m_d;
    parameters m_params;
    Eigen::SparseMatrix<real> m_weights;
    std::vector<real> diagonal;
    Eigen::SparseMatrix<real> m_diagonal;
    std::vector<real> m_eigenvectors;
    std::vector<size_t> m_groups;
    real m_duration;
  };

END_PROJECT_NAMESPACE
# endif
