#include <string>
#include <iostream>
#include <Eigen/Eigen>

using namespace std;

template <typename Scalar, typename Roots> inline void
	computeRoots (const Scalar& b, const Scalar& c, Roots& roots)
{
	//roots of quadratic equation
	roots (0) = Scalar (0);
	Scalar d = Scalar (b * b - 4.0 * c);
	if (d < 0.0)  // no real roots ! THIS SHOULD NOT HAPPEN!
		d = 0.0;

	Scalar sd = ::std::sqrt (d);

	roots (2) = 0.5f * (b + sd);
	roots (1) = 0.5f * (b - sd);
}

//////////////////////////////////////////////////////////////////////////////////////////
template <typename Matrix, typename Vector> inline void
	computeEigenvaluesforScaledMat (const Matrix& m, Vector& roots)
{
	typedef typename Matrix::Scalar Scalar;

	// The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
	// eigenvalues are the roots to this equation, all guaranteed to be
	// real-valued, because the matrix is symmetric.
	Scalar c0 =      m (0, 0) * m (1, 1) * m (2, 2)
		+ Scalar (2) * m (0, 1) * m (0, 2) * m (1, 2)
		- m (0, 0) * m (1, 2) * m (1, 2)
		- m (1, 1) * m (0, 2) * m (0, 2)
		- m (2, 2) * m (0, 1) * m (0, 1);
	Scalar c1 = m (0, 0) * m (1, 1) -
		m (0, 1) * m (0, 1) +
		m (0, 0) * m (2, 2) -
		m (0, 2) * m (0, 2) +
		m (1, 1) * m (2, 2) -
		m (1, 2) * m (1, 2);
	Scalar c2 = m (0, 0) + m (1, 1) + m (2, 2);

	if (fabs (c0) < Eigen::NumTraits < Scalar > ::epsilon ())  // one root is 0 -> quadratic equation
		computeRoots (c2, c1, roots);
	else
	{
		const Scalar s_inv3 = Scalar (1.0 / 3.0);
		const Scalar s_sqrt3 = std::sqrt (Scalar (3.0));
		// Construct the parameters used in classifying the roots of the equation
		// and in solving the equation for the roots in closed form.
		Scalar c2_over_3 = c2 * s_inv3;
		Scalar a_over_3 = (c1 - c2 * c2_over_3) * s_inv3;
		if (a_over_3 > Scalar (0))
			a_over_3 = Scalar (0);

		Scalar half_b = Scalar (0.5) * (c0 + c2_over_3 * (Scalar (2) * c2_over_3 * c2_over_3 - c1));

		Scalar q = half_b * half_b + a_over_3 * a_over_3 * a_over_3;
		if (q > Scalar (0))
			q = Scalar (0);

		// Compute the eigenvalues by solving for the roots of the polynomial.
		Scalar rho = std::sqrt (-a_over_3);
		Scalar theta = std::atan2 (std::sqrt (-q), half_b) * s_inv3;
		Scalar cos_theta = std::cos (theta);
		Scalar sin_theta = std::sin (theta);
		roots (0) = c2_over_3 + Scalar (2) * rho * cos_theta;
		roots (1) = c2_over_3 - rho * (cos_theta + s_sqrt3 * sin_theta);
		roots (2) = c2_over_3 - rho * (cos_theta - s_sqrt3 * sin_theta);

		// Sort in increasing order.
		if (roots (0) >= roots (1))
			std::swap (roots (0), roots (1));
		if (roots (1) >= roots (2))
		{
			std::swap (roots (1), roots (2));
			if (roots (0) >= roots (1))
				std::swap (roots (0), roots (1));
		}

		if (roots (0) <= 0)  // eigenval for symetric positive semi-definite matrix can not be negative! Set it to 0
			computeRoots (c2, c1, roots);
	}
}

template <typename Matrix, typename Vector> inline void
	computeEigenvalues (const Matrix& mat, Vector& evals)
{
	
	typedef typename Matrix::Scalar Scalar;
	// Scale the matrix so its entries are in [-1,1].  The scaling is applied
	// only when at least one matrix entry has magnitude larger than 1.

	Scalar scale = mat.cwiseAbs ().maxCoeff ();
	if (scale <= std::numeric_limits < Scalar > ::min ())
		scale = Scalar (1.0);

	Matrix scaledMat = mat / scale;
	// Compute the eigenvalues
	computeEigenvaluesforScaledMat (scaledMat, evals);
	evals *= scale;
}

template <typename Matrix, typename Vector> inline void
	computeEigen (const Matrix& mat, Matrix& evecs, Vector& evals)
{
	typedef typename Matrix::Scalar Scalar;
	// Scale the matrix so its entries are in [-1,1].  The scaling is applied
	// only when at least one matrix entry has magnitude larger than 1.

	Scalar scale = mat.cwiseAbs ().maxCoeff ();
	if (scale <= std::numeric_limits < Scalar > ::min ())
		scale = Scalar (1.0);

	Matrix scaledMat = mat / scale;

	// Compute the eigenvalues
	computeEigenvaluesforScaledMat (scaledMat, evals);

	if ( (evals (2) - evals (0)) <= Eigen::NumTraits < Scalar > ::epsilon ())
	{
		// all three equal
		evecs.setIdentity ();
	}
	else if ( (evals (1) - evals (0)) <= Eigen::NumTraits < Scalar > ::epsilon ())
	{
		// first and second equal
		Matrix tmp;
		tmp = scaledMat;
		tmp.diagonal ().array () -= evals (2);

		Vector vec1 = tmp.row (0).cross (tmp.row (1));
		Vector vec2 = tmp.row (0).cross (tmp.row (2));
		Vector vec3 = tmp.row (1).cross (tmp.row (2));

		Scalar len1 = vec1.squaredNorm ();
		Scalar len2 = vec2.squaredNorm ();
		Scalar len3 = vec3.squaredNorm ();

		if (len1 >= len2 && len1 >= len3)
			evecs.col (2) = vec1 / std::sqrt (len1);
		else if (len2 >= len1 && len2 >= len3)
			evecs.col (2) = vec2 / std::sqrt (len2);
		else
			evecs.col (2) = vec3 / std::sqrt (len3);

		evecs.col (1) = evecs.col (2).unitOrthogonal ();
		evecs.col (0) = evecs.col (1).cross (evecs.col (2));
	}
	else if ( (evals (2) - evals (1)) <= Eigen::NumTraits < Scalar > ::epsilon ())
	{
		// second and third equal
		Matrix tmp;
		tmp = scaledMat;
		tmp.diagonal ().array () -= evals (0);

		Vector vec1 = tmp.row (0).cross (tmp.row (1));
		Vector vec2 = tmp.row (0).cross (tmp.row (2));
		Vector vec3 = tmp.row (1).cross (tmp.row (2));

		Scalar len1 = vec1.squaredNorm ();
		Scalar len2 = vec2.squaredNorm ();
		Scalar len3 = vec3.squaredNorm ();

		if (len1 >= len2 && len1 >= len3)
			evecs.col (0) = vec1 / std::sqrt (len1);
		else if (len2 >= len1 && len2 >= len3)
			evecs.col (0) = vec2 / std::sqrt (len2);
		else
			evecs.col (0) = vec3 / std::sqrt (len3);

		evecs.col (1) = evecs.col (0).unitOrthogonal ();
		evecs.col (2) = evecs.col (0).cross (evecs.col (1));
	}
	else
	{
		Matrix tmp;
		tmp = scaledMat;
		tmp.diagonal ().array () -= evals (2);

		Vector vec1 = tmp.row (0).cross (tmp.row (1));
		Vector vec2 = tmp.row (0).cross (tmp.row (2));
		Vector vec3 = tmp.row (1).cross (tmp.row (2));

		Scalar len1 = vec1.squaredNorm ();
		Scalar len2 = vec2.squaredNorm ();
		Scalar len3 = vec3.squaredNorm ();
#ifdef _WIN32
		Scalar *mmax = new Scalar[3];
#else
		Scalar mmax[3];
#endif
		unsigned int min_el = 2;
		unsigned int max_el = 2;
		if (len1 >= len2 && len1 >= len3)
		{
			mmax[2] = len1;
			evecs.col (2) = vec1 / std::sqrt (len1);
		}
		else if (len2 >= len1 && len2 >= len3)
		{
			mmax[2] = len2;
			evecs.col (2) = vec2 / std::sqrt (len2);
		}
		else
		{
			mmax[2] = len3;
			evecs.col (2) = vec3 / std::sqrt (len3);
		}

		tmp = scaledMat;
		tmp.diagonal ().array () -= evals (1);

		vec1 = tmp.row (0).cross (tmp.row (1));
		vec2 = tmp.row (0).cross (tmp.row (2));
		vec3 = tmp.row (1).cross (tmp.row (2));

		len1 = vec1.squaredNorm ();
		len2 = vec2.squaredNorm ();
		len3 = vec3.squaredNorm ();
		if (len1 >= len2 && len1 >= len3)
		{
			mmax[1] = len1;
			evecs.col (1) = vec1 / std::sqrt (len1);
			min_el = len1 <= mmax[min_el] ? 1 : min_el;
			max_el = len1 > mmax[max_el] ? 1 : max_el;
		}
		else if (len2 >= len1 && len2 >= len3)
		{
			mmax[1] = len2;
			evecs.col (1) = vec2 / std::sqrt (len2);
			min_el = len2 <= mmax[min_el] ? 1 : min_el;
			max_el = len2 > mmax[max_el] ? 1 : max_el;
		}
		else
		{
			mmax[1] = len3;
			evecs.col (1) = vec3 / std::sqrt (len3);
			min_el = len3 <= mmax[min_el] ? 1 : min_el;
			max_el = len3 > mmax[max_el] ? 1 : max_el;
		}

		tmp = scaledMat;
		tmp.diagonal ().array () -= evals (0);

		vec1 = tmp.row (0).cross (tmp.row (1));
		vec2 = tmp.row (0).cross (tmp.row (2));
		vec3 = tmp.row (1).cross (tmp.row (2));

		len1 = vec1.squaredNorm ();
		len2 = vec2.squaredNorm ();
		len3 = vec3.squaredNorm ();
		if (len1 >= len2 && len1 >= len3)
		{
			mmax[0] = len1;
			evecs.col (0) = vec1 / std::sqrt (len1);
			min_el = len3 <= mmax[min_el] ? 0 : min_el;
			max_el = len3 > mmax[max_el] ? 0 : max_el;
		}
		else if (len2 >= len1 && len2 >= len3)
		{
			mmax[0] = len2;
			evecs.col (0) = vec2 / std::sqrt (len2);
			min_el = len3 <= mmax[min_el] ? 0 : min_el;
			max_el = len3 > mmax[max_el] ? 0 : max_el;
		}
		else
		{
			mmax[0] = len3;
			evecs.col (0) = vec3 / std::sqrt (len3);
			min_el = len3 <= mmax[min_el] ? 0 : min_el;
			max_el = len3 > mmax[max_el] ? 0 : max_el;
		}

		unsigned mid_el = 3 - min_el - max_el;
		evecs.col (min_el) = evecs.col ( (min_el + 1) % 3).cross (evecs.col ( (min_el + 2) % 3)).normalized ();
		evecs.col (mid_el) = evecs.col ( (mid_el + 1) % 3).cross (evecs.col ( (mid_el + 2) % 3)).normalized ();
#ifdef _WIN32
		delete [] mmax;
#endif
	}
	// Rescale back to the original size.
	evals *= scale;
}
template <typename PointT,typename Matrix> inline unsigned int
	filter (const pcl::PointCloud<PointT> &cloud,
	const std::vector<int> &indices,	
	Matrix &result
	)
{	
	result.setZero(indices.size(),3);
	int count=0;
	if (cloud.is_dense)
	{
		for (std::vector<int>::const_iterator it = indices.begin (); it != indices.end (); ++it)
		{
			result(count,0)=cloud[*it].x;
			result(count,1)=cloud[*it].y;
			result(count,2)=cloud[*it].z;
			++count;
		}
	}
	else
	{
		for (std::vector<int>::const_iterator it = indices.begin (); it != indices.end (); ++it)
		{
			if (!isFinite (cloud[*it]))
				continue;
			result(count,0)=cloud[*it].x;
			result(count,1)=cloud[*it].y;
			result(count,2)=cloud[*it].z;
			++count;
		}
		result=result.topRows(count).eval();
	}
	return count;
}

template <typename MatrixIn, typename MatrixOut> inline void
	computeCovariance (const MatrixIn& mat, MatrixOut& covMat)
{	
	// compute Covariance using Eigen
	covMat.setZero ();
	//at least 3 points are required
	if (mat.rows()>=3)	
	{
		typedef typename MatrixIn::Scalar Scalar;	
		MatrixIn input = mat;
		Eigen::Matrix< Scalar,1,Eigen::Dynamic> meanVal=input.colwise().mean(); 
		input.rowwise() -= meanVal;	
		covMat= input.adjoint() * input /(input.rows()-1);  
	}
	return;
}
template <typename MatrixIn, typename MatrixOut,typename VectorOut> inline void
	computeCovariance (
	const MatrixIn & mat, 
	MatrixOut & covMat,
	VectorOut & xyz_centroid)
{			
	// compute Covariance and xyz_centroid using Eigen
	covMat.setZero ();
	xyz_centroid.setZero ();
	//at least 3 points are required
	if (mat.rows()>=3)	
	{
		typedef typename MatrixIn::Scalar Scalar;
		MatrixIn input = mat;
		Eigen::Matrix< Scalar,1,Eigen::Dynamic> meanVal=input.colwise().mean(); 
		input.rowwise() -= meanVal;	
		covMat= input.adjoint() * input /(input.rows()-1);  
		xyz_centroid=VectorOut::Map(meanVal.data());
	}
	return;	
}
