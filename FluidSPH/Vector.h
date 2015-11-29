/* $Id: Vector.h,v 1.1.1.1 2005/05/01 15:18:55 ovidiom Exp $ */

#ifndef _VECTOR2D_H_
#define _VECTOR2D_H_

#include <cmath>


struct Vector2d
{
	double	x, y;

	inline Vector2d(){
		x = y = 0.0;
	}
	
	inline Vector2d(double x, double y){
		this->x = x;
		this->y = y;
	}
	
	inline Vector2d(double xy){
		x = y = xy;
	}

	inline Vector2d(const double *xyArr){
		x = xyArr[0];
		y = xyArr[1];
	}
	
	inline void operator +=(double s){
		x += s;
		y += s;
	}
	
	inline void operator +=(const Vector2d &v){
		x += v.x;
		y += v.y;
	}
	
	inline void operator -=(double s){
		x -= s;
		y -= s;
	}
	
	inline void operator -=(const Vector2d &v){
		x -= v.x;
		y -= v.y;
	}
	
	inline void operator *=(double s){
		x *= s;
		y *= s;
	}
	
	inline void operator *=(const Vector2d &v){
		x *= v.x;
		y *= v.y;
	}
	
	inline void operator /=(double s){
		double	inv = 1.0 / s;

		x *= inv;
		y *= inv;
	}

	inline void operator /=(const Vector2d &v){
		x /= v.x;
		y /= v.y;
	}
};


inline Vector2d operator +(const Vector2d &v, double s){
	return (Vector2d(v.x + s, v.y + s));
}

inline Vector2d operator +(double s, const Vector2d &v){
	return (Vector2d(s + v.x, s + v.y));
}

inline Vector2d operator +(const Vector2d &u, const Vector2d &v){
	return (Vector2d(u.x + v.x, u.y + v.y));
}

inline Vector2d operator -(const Vector2d &v, double s){
	return (Vector2d(v.x - s, v.y - s));
}

inline Vector2d operator -(double s, const Vector2d &v){
	return (Vector2d(s - v.x, s - v.y));
}

inline Vector2d operator -(const Vector2d &u, const Vector2d &v){
	return (Vector2d(u.x - v.x, u.y - v.y));
}

inline Vector2d operator *(const Vector2d &v, double s){
	return (Vector2d(v.x * s, v.y * s));
}

inline Vector2d operator *(double s, const Vector2d &v){
	return (Vector2d(s * v.x, s * v.y));
}

inline Vector2d operator *(const Vector2d &u, const Vector2d &v){
	return (Vector2d(u.x * v.x, u.y * v.y));
}

inline Vector2d operator /(const Vector2d &v, double s){
	double	inv = 1.0f / s;

	return (Vector2d(v.x * inv, v.y * inv));
}

inline Vector2d operator /(double s, const Vector2d &v){
	return (Vector2d(s / v.x, s / v.y));
}

inline Vector2d operator /(const Vector2d &u, const Vector2d &v){
	return (Vector2d(u.x / v.x, u.y / v.y));
}

inline Vector2d operator -(const Vector2d &v){
	return (Vector2d(-v.x, -v.y));
}

inline double dot(const Vector2d &u, const Vector2d &v){
	return (u.x * v.x + u.y * v.y);
}

inline double length(const Vector2d &v)
{
	return sqrt(v.x * v.x + v.y * v.y);
}

inline Vector2d normalize(const Vector2d &v){
	return (v / length(v));
}



#endif /* _VECTOR2D_H_ */

