#ifndef MY_RAND_HPP
#define MY_RAND_HPP

#include <random>

class my_random{

public:

	static std::default_random_engine &get_gre(int seed = 0)
	{
		static std::default_random_engine gre(seed);
		return gre;
	}
	
	static void set_seed(int s)
	{
		get_gre().seed(s);
	}
	
private:
	my_random(){}
	my_random(const my_random &);
	void operator=(const my_random &);
};



template <typename real_type> inline
real_type randf(real_type lo = 0.0, real_type hi = 1.0)
{
	std::uniform_real_distribution<real_type> urd(lo,hi);
	return urd( my_random::get_gre() );
}

inline
double randf(double lo = 0.0, double hi = 1.0)
{
  return randf<double>(lo,hi);
}

// This is NON-inclusive!
template <typename int_type> inline
int_type randi( int_type lo, int_type hi )
{
	std::uniform_int_distribution<int_type> urd(lo,hi-1);
	return urd( my_random::get_gre() );
}


template <typename real_type> inline
real_type randn(real_type mu=0.0, real_type sigma=1.0)
{
	std::normal_distribution<real_type> nd(mu,sigma);
	return nd( my_random::get_gre() );
}



inline bool flip_coin()
{
	std::uniform_int_distribution<int> urd(0,1);
	return urd( my_random::get_gre() );
}



#endif // MY_RAND_HPP