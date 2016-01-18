/*  Created on: Dec 29, 2015
 *      Author: T. Delame (tdelame@gmail.com)
 */

# ifndef OMP_KMEAN_H_
# define OMP_KMEAN_H_
# include <project.h>
# include <vector>
# include <functional>
BEGIN_PROJECT_NAMESPACE 

  class omp_kmean {
  public:
    omp_kmean(
        real* input,
        size_t number_of_points,
        size_t dimension_of_points,
        size_t number_of_clusters,
        std::vector<size_t>& groups);

    void execute();

    real
    get_execution_time() const;

    const real*
    get_centroid( size_t cluster_id ) const;

    size_t
    get_group( size_t point_id ) const;

    size_t
    get_number_of_clusters() const;

    size_t
    get_point_dimension() const;

    size_t
    get_number_of_points() const;

  protected:
    void init_centroids();
    void recompute_centroids();
    void recompute_clusters();

    typedef std::function<real(real*,real*)> square_distance_function;
    square_distance_function m_sdistance;

    std::vector<size_t>& m_groups;
    std::vector<real> m_centroids;
    real* m_input;
    real m_duration;
    const size_t m_n;
    const size_t m_d;
    const size_t m_k;
    size_t m_points_changed;
  };

 END_PROJECT_NAMESPACE

# endif
