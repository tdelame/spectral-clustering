/*  Created on: Dec 29, 2015
 *      Author: T. Delame (tdelame@gmail.com)
 */
# include <omp_kmean.h>
# include <algorithm>
# include <omp.h>
BEGIN_PROJECT_NAMESPACE

  omp_kmean::omp_kmean(
      real* input, size_t number_of_points,
      size_t dimension_of_points, size_t number_of_clusters,
      std::vector<size_t>& groups)
  : m_groups{ groups }, m_centroids (number_of_clusters * dimension_of_points, 0),
    m_input{ input }, m_duration{ 0 }, m_n{ number_of_points }, m_d{ dimension_of_points },
    m_k{ number_of_clusters }, m_points_changed{ 0 }
  {
    if (!m_d)
      {
        LOG(fatal, "the dimension of points cannot be 0");
        exit ( EXIT_FAILURE);
      }
    else if (m_d == 1)
      {
        m_sdistance = []( real* a, real* b ) -> real
          {
            real diff = *a - *b;
            return diff * diff;
          };
      }
    else if (m_d == 2)
      {
        m_sdistance = []( real* a, real* b ) -> real
          {
            real diff0 = a[0] - b[0];
            real diff1 = a[1] - b[1];
            return diff0 * diff0 + diff1 * diff1;
          };
      }
    else if (m_d == 3)
      {
        m_sdistance = []( real* a, real* b ) -> real
          {
            real diff0 = a[0] - b[0];
            real diff1 = a[1] - b[1];
            real diff2 = a[2] - b[2];
            return diff0 * diff0 + diff1 * diff1 + diff2 * diff2;
          };
      }
    else
      {
        // this is adapted from the flann::L2<real> structure
        m_sdistance = [dimension_of_points]( real* a, real* b) -> real
          {
            real result = 0;
            real diff0, diff1, diff2, diff3;
            real* last = a + dimension_of_points;
            real* lastgroup = last - 3;

            /* Process 4 items with each loop for efficiency. */
            while (a < lastgroup)
              {
                diff0 = real(a[0] - b[0]);
                diff1 = real(a[1] - b[1]);
                diff2 = real(a[2] - b[2]);
                diff3 = real(a[3] - b[3]);
                result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
                a += 4;
                b += 4;
              }
            /* Process last 0-3 pixels.  Not needed for standard vector lengths. */
            while (a < last)
              {
                diff0 = real(*a++ - *b++);
                result += diff0 * diff0;
              }
            return result;
          };
      }
    m_groups.resize( number_of_points);
    std::fill( m_groups.begin(), m_groups.end(), m_k );
  }

  void
  omp_kmean::execute ()
  {
    m_duration = omp_get_wtime ();
    init_centroids ();
    do
      {
        recompute_centroids ();
        recompute_clusters ();
        std::cout << "one k-mean iteration"<<std::endl;
      }
    // stops when 99.9% of the points do not change clusters
    // (m_n >> 10 equals m_n / 1024, thus represents roughly 0.1% of m_n)
    while ( m_points_changed > (m_n >> 10));
    m_duration = omp_get_wtime () - m_duration;
  }

  real
  omp_kmean::get_execution_time() const
  {
    return m_duration;
  }

  const real*
  omp_kmean::get_centroid (size_t cluster_id) const
  {
    return &m_centroids[cluster_id * m_d];
  }

  size_t
  omp_kmean::get_group (size_t point_id) const
  {
    return m_groups[point_id];
  }

  size_t
  omp_kmean::get_number_of_clusters () const
  {
    return m_k;
  }

  size_t
  omp_kmean::get_point_dimension () const
  {
    return m_d;
  }

  void
  omp_kmean::init_centroids ()
  {
    // the first centroid is located at a random input point
      {
        real* source = m_input + size_t (unit_random () * m_n) * m_d;
        // centroid_0 = source
        std::copy (source, source + m_d, m_centroids.data ());
      }

    // compute the other centroids with a kmean ++ algorithm
    std::vector<real> sdistances (m_n, real (0));
    for (size_t k = 1; k < m_k; ++k)
      {
        real sum = 0;
        # pragma omp parallel
          {
            real thread_sum = 0;
            # pragma omp for schedule(static)
            for (size_t i = 0; i < m_n; ++i)
              {
                real* point_i = m_input + m_d * i;
                real* centroid_j = m_centroids.data ();
                auto min_sdi = m_sdistance (centroid_j, point_i);
                for (size_t j = 1; j < k; ++j, centroid_j += m_d)
                  {
                    min_sdi = std::min (min_sdi,
                                        m_sdistance (centroid_j, point_i));
                  }
                sdistances[i] = min_sdi;
                thread_sum += min_sdi;
              }
            # pragma omp critical
            sum += thread_sum;
          }
        sum *= unit_random ();

        auto point_i = m_input;
        for (size_t i = 0; i < m_n; ++i, point_i += m_d)
          {
            if ((sum -= sdistances[i]) > 0)
              continue;
            // centroid_k = point_i
            std::copy (point_i, point_i + m_d, m_centroids.data () + k * m_d);
            break;
          }
      }
    recompute_clusters ();
  }

  void
  omp_kmean::recompute_centroids ()
  {
    std::vector<size_t> counts (m_k, 0);
    // centroid_i = 0
    std::fill (m_centroids.begin (), m_centroids.end (), real (0));

// not applied as it make the run 6 time slower! The reason might be the cache use.
//    # pragma omp parallel for
//    for( size_t i = 0; i < m_n; ++ i )
//      {
//        auto g = m_groups[ i ];
//        auto centroid = m_centroids.data() + m_d * g;
//        real* point_i = m_input + m_d * i;
//        std::transform( point_i, point_i + m_d, centroid, centroid, std::plus<real>());
//        # pragma omp critical
//        ++counts[ g ];
//      }

    real* point_i = m_input;
    for (size_t i = 0; i < m_n; ++i)
      {
        auto g = m_groups[i];
        auto centroid = m_centroids.data () + m_d * g;

        for (size_t d = 0; d < m_d; ++d, ++centroid, ++point_i)
          *centroid += *point_i;
      }

    for (auto& g : m_groups)
      ++counts[g];

    // no impact at all
//    # pragma omp parallel for
    for (size_t k = 0; k < m_k; ++k)
      {
        const auto divisor = real (1.0) / real (counts[k]);
        auto centroid_k = m_centroids.data () + m_d * k;
        // centroid_k *= divisor
        for (size_t d = 0; d < m_d; ++d, ++centroid_k)
          *centroid_k *= divisor;
      }
  }

  void
  omp_kmean::recompute_clusters ()
  {
    m_points_changed = 0;

    # pragma omp parallel for
    for (size_t i = 0; i < m_n; ++i)
      {
        size_t g = 0;
        real mind2 = REAL_MAX;
        auto centroid_k = m_centroids.data ();
        auto point_i = m_input + i * m_d;
        for (size_t k = 0; k < m_k; ++k, centroid_k += m_d)
          {
            auto d2 = m_sdistance (point_i, centroid_k);
            if (d2 < mind2)
              {
                mind2 = d2;
                g = k;
              }
          }

        if (g != m_groups[i])
          {
            # pragma omp critical
            ++m_points_changed;
          }
        m_groups[i] = g;
      }
  }

  size_t
  omp_kmean::get_number_of_points() const
  {
    return m_n;
  }

END_PROJECT_NAMESPACE
