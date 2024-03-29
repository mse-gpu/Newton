struct vector{
	float x;
	float y;
};

#define LIMIT 1000
#define PRECISION 1
#define CIRCLE 3
#define SQRT3 1.7320508075688772935

bool near(float src, float target){
    float delta = src - target;

    if(delta < 0){
	delta = -delta;
    }

    if(delta <= PRECISION){
	return true;
    }

    return false;
}

inline int real_newton(float x, float y){
    vector xn = {x,y};

    int current = 0;

    int times = 0;
    int last = 0;

    while(current < LIMIT){
	float fnx = xn.x * xn.x * xn.x - 3 * xn.x * xn.y * xn.y - 1;
	float fny = xn.y * xn.y * xn.y - 3 * xn.x * xn.x * xn.y;

	float ja = 3 * xn.x * xn.x - 3 * xn.y * xn.y;
	float jd = 3 * xn.y * xn.y - 3 * xn.x * xn.x;
	float jbc = 6 * xn.x * xn.y;

	float det = ja * jd - jbc * jbc; //det(A) = a*d - b*c

	float dx = (jd / det) * fnx + (jbc / det) * fny;
	float dy = (jbc / det) * fnx + (ja / det) * fny;

	xn.x = xn.x - dx;
	xn.y = xn.y - dy;

	if(near(xn.x, 1) && near(xn.y, 0)){
	    if(times == CIRCLE && last == 1){
		return 1;
	    }

	    if(last == 1){
		++times;
	    } else {
		times = 1;
	    }

	    last = 1;
	} else if(near(xn.x, -1/2) && near(xn.y, SQRT3 / 2)){
	    if(times == CIRCLE && last == 2){
		return 2;
	    }

	    if(last == 2){
		++times;
	    } else {
		times = 1;
	    }

	    last = 2;
	} else if(near(xn.x, -1/2) && near(xn.y, -SQRT3 / 2)){
	    if(times == CIRCLE && last == 3){
		return 3;
	    }

	    if(last == 3){
		++times;
	    } else {
		times = 1;
	    }

	    last = 3;
	} else {
	    times = 0;
	    last = 0;
	}

	++current;
    }

    //Once we are here, it means that we are out the loop: black point
    return 0;
}
