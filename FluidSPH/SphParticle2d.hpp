#ifndef _POINT3F_HPP_
#define _POINT3F_HPP_
namespace sph
{
	class SphParticle2d
	{
	public:
		double x,y;
		double vx,vy;
		double ax,ay;
		double cr,cg,cb;
		double density;
		double pressure;
		SphParticle2d()
		{
			vy=vx=0;
			ax=ay=0;
			cr=cg=cb;
			density=pressure=0;
		}
		~SphParticle2d()
		{
		}
	};
}
#endif