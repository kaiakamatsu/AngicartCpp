

// bitmaps for representing strings of numbers
GLubyte bitZero16[] = {0, 0, 3, 192, 15, 240, 28, 56, 48, 12, 48, 12, 96, 6, 96, 6,
		96, 6, 96, 6, 48, 12, 48, 12, 28, 56, 15, 240, 3, 192, 0, 0},
bitOne16[] = {0, 0, 0, 192, 0, 192, 0, 192, 0, 192, 0, 192, 0, 192, 0, 192,
		0, 192, 0, 192, 0, 192, 0, 192, 3, 192, 1, 192, 0, 192, 0, 0},
bitTwo16[] = {0, 0, 127, 254, 127, 254, 224, 0, 60, 0, 31, 0, 7, 224, 1, 248,
		0, 60, 0, 28, 96, 6, 96, 6, 56, 28, 31, 248, 7, 224, 0, 0},
bitThree16[] = {0, 0, 15, 248, 31, 252, 48, 14, 96, 6, 96, 6, 0, 28, 1, 248,
		1, 248, 0, 28, 96, 6, 96, 6, 48, 14, 31, 252, 15, 248, 0, 0},
bitFour16[] = {0, 0, 0, 24, 0, 24, 0, 24, 0, 24, 0, 24, 0, 24, 127, 254,
		127, 254, 96, 24, 96, 24, 96, 24, 96, 24, 96, 24, 96, 24, 0, 0},
bitFive16[] = {0, 0, 31, 248, 63, 252, 48, 14, 96, 6, 96, 6, 0, 6, 0, 6,
		112, 62, 127, 248, 111, 224, 96, 0, 96, 0, 127, 254, 127, 254, 0, 0},
bitSix16[] = {0, 0, 3, 240, 15, 248, 28, 28, 60, 6, 108, 6, 110, 6, 103, 252,
		99, 248, 96, 0, 48, 0, 48, 0, 28, 0, 15, 192, 3, 192, 0, 0},
bitSeven16[] = {0, 0, 1, 192, 1, 192, 0, 192, 0, 224, 0, 96, 0, 112, 0, 112, 0, 48,
		0, 56, 0, 24, 0, 28, 0, 12, 0, 14, 127, 254, 127, 254, 0, 0},
bitEight16[] = {0, 0, 7, 224, 31, 248, 56, 28, 96, 6, 96, 6, 56, 28, 31, 248,
		31, 248, 56, 28, 96, 6, 96, 6, 56, 28, 31, 248, 7, 224, 0, 0},
bitNine16[] = {0, 0, 0, 6, 0, 6, 0, 6, 0, 6, 0, 6, 0, 6, 31, 254,
		63, 254, 112, 14, 96, 6, 96, 6, 112, 14, 63, 254, 31, 252, 0, 0};
GLubyte bitZero8[] = {0, 60, 102, 66, 66, 102, 60, 0},
bitOne8[] = {0, 8, 8, 8, 8, 24, 8, 0},
bitTwo8[] = {0, 126, 32, 28, 2, 66, 60, 0},
bitThree8[] = { 0, 60, 66, 2, 12, 66, 60, 0},
bitFour8[] = { 0, 4, 4, 4, 126, 68, 68, 0},
bitFive8[] = { 0, 124, 2, 2, 124, 64, 126, 0},
bitSix8[] = { 0, 60, 114, 94, 64, 64, 56, 0},
bitSeven8[] = { 0, 16, 16, 8, 4, 2, 126, 0},
bitEight8[] = {0, 60, 66, 66, 60, 66, 60, 0},
bitNine8[] = {0, 2, 2, 62, 66, 66, 60, 0};

GLubyte bitSpace16f[64], bitPeriod16f[64], bitDash16f[64], bitZero16f[64], bitOne16f[64], bitTwo16f[64], bitThree16f[64], bitFour16f[64], bitFive16f[64], bitSix16f[64], bitSeven16f[64], bitEight16f[64], bitNine16f[64];
GLubyte bitSpace8f[32], bitPeriod8f[32], bitDash8f[32], bitZero8f[32], bitOne8f[32], bitTwo8f[32], bitThree8f[32], bitFour8f[32], bitFive8f[32], bitSix8f[32], bitSeven8f[32], bitEight8f[32], bitNine8f[32];

void initBitmapNumbers(){
	for(unsigned int i(0); i < 16; i++){
		for(unsigned int j(0); j < 2; j++){
			bitZero16f[4*i + j] = bitZero16[2*i + j];
			bitOne16f[4*i + j] = bitOne16[2*i + j];
			bitTwo16f[4*i + j] = bitTwo16[2*i + j];
			bitThree16f[4*i + j] = bitThree16[2*i + j];
			bitFour16f[4*i + j] = bitFour16[2*i + j];
			bitFive16f[4*i + j] = bitFive16[2*i + j];
			bitSix16f[4*i + j] = bitSix16[2*i + j];
			bitSeven16f[4*i + j] = bitSeven16[2*i + j];
			bitEight16f[4*i + j] = bitEight16[2*i + j];
			bitNine16f[4*i + j] = bitNine16[2*i + j];
		}
		bitZero16f[4*i + 2] = bitZero16f[4*i + 3] = 0;
		bitOne16f[4*i + 2] = bitOne16f[4*i + 3] = 0;
		bitTwo16f[4*i + 2] = bitTwo16f[4*i + 3] = 0;
		bitThree16f[4*i + 2] = bitThree16f[4*i + 3] = 0;
		bitFour16f[4*i + 2] = bitFour16f[4*i + 3] = 0;
		bitFive16f[4*i + 2] = bitFive16f[4*i + 3] = 0;
		bitSix16f[4*i + 2] = bitSix16f[4*i + 3] = 0;
		bitSeven16f[4*i + 2] = bitSeven16f[4*i + 3] = 0;
		bitEight16f[4*i + 2] = bitEight16f[4*i + 3] = 0;
		bitNine16f[4*i + 2] = bitNine16f[4*i + 3] = 0;
	}
	for(unsigned int i(0); i < 64; i++)
		bitSpace16f[i] = bitPeriod16f[i] = bitDash16f[i] = 0;
	bitPeriod16f[4] = bitPeriod16f[8] = bitPeriod16f[12]= 3;
	bitPeriod16f[5] = bitPeriod16f[9] = bitPeriod16f[13] = 192;
	bitDash16f[28] = bitDash16f[32] = 63;
	bitDash16f[29] = bitDash16f[33] = 252;
	
	for(unsigned int i(0); i < 8; i++){
		bitZero8f[4*i] = bitZero8[i];
		bitOne8f[4*i] = bitOne8[i];
		bitTwo8f[4*i] = bitTwo8[i];
		bitThree8f[4*i] = bitThree8[i];
		bitFour8f[4*i] = bitFour8[i];
		bitFive8f[4*i] = bitFive8[i];
		bitSix8f[4*i] = bitSix8[i];
		bitSeven8f[4*i] = bitSeven8[i];
		bitEight8f[4*i] = bitEight8[i];
		bitNine8f[4*i] = bitNine8[i];
		for(unsigned int j(1); j < 4; j++){
			bitZero8f[4*i + j] = 0;
			bitOne8f[4*i + j] = 0;
			bitTwo8f[4*i + j] = 0;
			bitThree8f[4*i + j] = 0;
			bitFour8f[4*i + j] = 0;
			bitFive8f[4*i + j] = 0;
			bitSix8f[4*i + j] = 0;
			bitSeven8f[4*i + j] = 0;
			bitEight8f[4*i + j] = 0;
			bitNine8f[4*i + j] = 0;
		}
	}
	for(unsigned int i(0); i < 16; i++)
		bitSpace8f[i/2] = bitPeriod8f[i] = bitDash8f[i] = 0;
	bitPeriod8f[2] = bitPeriod8f[4] = 24;
	bitDash8f[10] = 60;
}

