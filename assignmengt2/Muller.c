float rtmuller(float (*func)(float), float x1, float x2, float xacc){
    float x3 = 0.5 * (x1 + x2);
    float x4;
    float a, b, c;
	const double eps = 1.084202e-19;
    do{
        c = func(x3);
		b = (x1 - x3) / (x2 - x3 + eps) / (x1 - x2 + eps) * (func(x2) - func(x3)) - (x2 - x3) / (x1 - x3 + eps) / (x1 - x2 + eps) * (func(x1) - func(x3));
        a = (func(x1) - func(x3)) / (x1 - x3 + eps) / (x1 - x2 + eps) - (func(x2) - func(x3)) / (x2 - x3 + eps) / (x1 - x2 + eps);
		if(b * b - 4 * a * c < 0){
			break;
		}
        x4 = x3 - 2 * c / (b + (b > 0 ? 1 : -1) * sqrt(b * b - 4 * a * c) + eps);
        x1 = x2; x2 = x3; x3 = x4;
    }while(fabs(func(x4)) > xacc);
    return x4;
}