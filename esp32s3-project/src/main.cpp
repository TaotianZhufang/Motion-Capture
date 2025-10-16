#include <SPI.h>
#include <Arduino.h>

//BMI1 FSPI
#define BMI1_MISO 9
#define BMI1_MOSI 11
#define BMI1_SCK  12
#define BMI1_CS   13

//BMI2 HSPI
#define BMI2_MISO 14
#define BMI2_MOSI 17
#define BMI2_SCK  18
#define BMI2_CS   21

SPIClass spi1(FSPI); //BMI1
SPIClass spi2(HSPI); //BMI2

//BMI323 Register
#define CHIP_ID_REG 0x00
#define ACC_CONF_REG 0x20
#define GYR_CONF_REG 0x21
#define ACC_DATA_X_REG 0x03
#define ACC_DATA_Y_REG 0x04
#define ACC_DATA_Z_REG 0x05
#define GYR_DATA_X_REG 0x06
#define GYR_DATA_Y_REG 0x07
#define GYR_DATA_Z_REG 0x08
#define TEMP_DATA_REG 0x09
#define CMD_REG 0x7E
#define FEATURE_IO0_REG 0x10
#define FEATURE_IO1_REG 0x11
#define FEATURE_IO2_REG 0x12
#define FEATURE_IO_STATUS_REG 0x14
#define FEATURE_CTRL_REG 0x40

//常量配置
#define SPI_FREQ 1000000
#define OUTPUT_RATE_HZ 20
#define ACC_RANGE_LSB_PER_16G 2048.0
#define GYR_RANGE_LSB_PER_2000DPS 16.384
#define SOFT_RESET_CMD 0xDEAF
#define ACC_CONF_NORMAL_100HZ_16G 0x4038
#define GYR_CONF_NORMAL_100HZ_2000DPS 0x4048

//Settings For BMI323 More In Bmi323Config.h
#define ACC_RANGE_LSB_PER_G ACC_RANGE_LSB_PER_16G
#define GYR_RANGE_LSB_PER_DPS GYR_RANGE_LSB_PER_2000DPS
#define ACC_CONF_NORMAL_HZ_G ACC_CONF_NORMAL_100HZ_16G
#define GYR_CONF_NORMAL_HZ_DPS GYR_CONF_NORMAL_100HZ_2000DPS

// System constants
#define deltat 0.05f                                    // sampling period in seconds (shown as 1 ms)
#define gyroMeasError 3.14159265358979 * (5000.0f / 180.0f) // gyroscope measurement error in rad/s (shown as 5 deg/s)
#define gyroMeasDrift 3.14159265358979 * (2000.0f / 180.0f) // gyroscope measurement error in rad/s/s (shown as 0.2f deg/s/s)
//#define beta sqrt(3.0f / 4.0f) * gyroMeasError           // compute beta
//#define zeta sqrt(3.0f / 4.0f) * gyroMeasDrift           // compute zeta
#define beta         0.03        // compute beta  较大的beta增强加速度计的修正，但可能引入运动干扰；
                                               // 较小的beta更依赖陀螺仪，但可能累积漂移。
#define zeta         0.0001      // compute zeta
#define Gyro_Gr	    0.01745329f				//角速度变成弧度(3.1415/180)
#define G			      9.7970f		        	//m/s^2
#define RadtoDeg    57.324841f				//弧度到角度 (弧度 * 180/3.1415)
#define DegtoRad    0.0174533f				//角度到弧度 (角度 * 3.1415/180)

// Global system variables
float a_x, a_y, a_z;
// accelerometer measurements
float w_x, w_y, w_z;
// gyroscope measurements in rad/s
float m_x=0, m_y=0, m_z=0;
// magnetometer measurements
volatile float SEq_1 = 1.00f, SEq_2 = 0.00f, SEq_3 = 0.00f, SEq_4 = 0.00f;
// estimated orientation quaternion elements with initial conditions
float b_x = 1, b_z = 0;
// reference direction of flux in earth frame
float w_bx = 0, w_by = 0, w_bz = 0;
// estimate gyroscope biases error

// Function to compute one filter iteration
void filterUpdate(float w_x, float w_y, float w_z, float a_x, float a_y, float a_z, float m_x, float m_y, float m_z)
{
// local system variables
    float norm;
// vector norm
    float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3, SEqDot_omega_4; // quaternion rate from gyroscopes elements
    float f_1, f_2, f_3, f_4, f_5, f_6;
// objective function elements
    float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33,
// objective function Jacobian elements
    J_41, J_42, J_43, J_44, J_51, J_52, J_53, J_54, J_61, J_62, J_63, J_64; //
    float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4;
// estimated direction of the gyroscope error
    float w_err_x, w_err_y, w_err_z;
// estimated direction of the gyroscope error (angular)
    float h_x, h_y, h_z;
// computed flux in the earth frame

// axulirary variables to avoid reapeated calcualtions
    float halfSEq_1 = 0.5f * SEq_1;
    float halfSEq_2 = 0.5f * SEq_2;
    float halfSEq_3 = 0.5f * SEq_3;
    float halfSEq_4 = 0.5f * SEq_4;

    float twoSEq_1 = 2.0f * SEq_1;
    float twoSEq_2 = 2.0f * SEq_2;
    float twoSEq_3 = 2.0f * SEq_3;
    float twoSEq_4 = 2.0f * SEq_4;

    float twob_x = 2.0f * b_x;
    float twob_z = 2.0f * b_z;

    float twob_xSEq_1 = 2.0f * b_x * SEq_1;
    float twob_xSEq_2 = 2.0f * b_x * SEq_2;
    float twob_xSEq_3 = 2.0f * b_x * SEq_3;
    float twob_xSEq_4 = 2.0f * b_x * SEq_4;
    float twob_zSEq_1 = 2.0f * b_z * SEq_1;
    float twob_zSEq_2 = 2.0f * b_z * SEq_2;
    float twob_zSEq_3 = 2.0f * b_z * SEq_3;
    float twob_zSEq_4 = 2.0f * b_z * SEq_4;

    float SEq_1SEq_2;
    float SEq_1SEq_3 = SEq_1 * SEq_3;
    float SEq_1SEq_4;
    float SEq_2SEq_3;
    float SEq_2SEq_4 = SEq_2 * SEq_4;
    float SEq_3SEq_4;
    float twom_x = 2.0f * m_x;
    float twom_y = 2.0f * m_y;
    float twom_z = 2.0f * m_z;

// normalise the accelerometer measurement
    norm = sqrt(a_x * a_x + a_y * a_y + a_z * a_z);
    a_x /= norm;
    a_y /= norm;
    a_z /= norm;
// normalise the magnetometer measurement
 //   norm = sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
  //  m_x /= norm;
  //  m_y /= norm;
  //  m_z /= norm;

// compute the objective function and Jacobian
    f_1 = twoSEq_2 * SEq_4 - twoSEq_1 * SEq_3 - a_x;
    f_2 = twoSEq_1 * SEq_2 + twoSEq_3 * SEq_4 - a_y;
    f_3 = 1.0f - twoSEq_2 * SEq_2 - twoSEq_3 * SEq_3 - a_z;
    //加速度计误差矩阵4*3

    f_4 = twob_x * (0.5f - SEq_3 * SEq_3 - SEq_4 * SEq_4) + twob_z * (SEq_2SEq_4 - SEq_1SEq_3) - m_x;
    f_5 = twob_x * (SEq_2 * SEq_3 - SEq_1 * SEq_4) + twob_z * (SEq_1 * SEq_2 + SEq_3 * SEq_4) - m_y;
    f_6 = twob_x * (SEq_1SEq_3 + SEq_2SEq_4) + twob_z * (0.5f - SEq_2 * SEq_2 - SEq_3 * SEq_3) - m_z;
    //磁力计误差矩阵4*3

    J_11or24 = twoSEq_3;
// J_11 negated in matrix multiplication
    J_12or23 = 2.0f * SEq_4;
    J_13or22 = twoSEq_1;
// J_12 negated in matrix multiplication
    J_14or21 = twoSEq_2;
    J_32 = 2.0f * J_14or21;
// negated in matrix multiplication
    J_33 = 2.0f * J_11or24;
// negated in matrix multiplication
    J_41 = twob_zSEq_3;
// negated in matrix multiplication
    J_42 = twob_zSEq_4;
    J_43 = 2.0f * twob_xSEq_3 + twob_zSEq_1;
// negated in matrix multiplication
    J_44 = 2.0f * twob_xSEq_4 - twob_zSEq_2;
// negated in matrix multiplication
    J_51 = twob_xSEq_4 - twob_zSEq_2;
// negated in matrix multiplication
    J_52 = twob_xSEq_3 + twob_zSEq_1;
    J_53 = twob_xSEq_2 + twob_zSEq_4;
    J_54 = twob_xSEq_1 - twob_zSEq_3;
// negated in matrix multiplication
    J_61 = twob_xSEq_3;
    J_62 = twob_xSEq_4 - 2.0f * twob_zSEq_2;
    J_63 = twob_xSEq_1 - 2.0f * twob_zSEq_3;
    J_64 = twob_xSEq_2;
//6*4的Jacobian矩阵的每一式

// compute the gradient (matrix multiplication)
    SEqHatDot_1 = J_14or21 * f_2 - J_11or24 * f_1 - J_41 * f_4 - J_51 * f_5 + J_61 * f_6;
    SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2 - J_32 * f_3 + J_42 * f_4 + J_52 * f_5 + J_62 * f_6;
    SEqHatDot_3 = J_12or23 * f_2 - J_33 * f_3 - J_13or22 * f_1 - J_43 * f_4 + J_53 * f_5 + J_63 * f_6;
    SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2 - J_44 * f_4 - J_54 * f_5 + J_64 * f_6;

// normalise the gradient to estimate direction of the gyroscope error
    norm = sqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3 + SEqHatDot_4 * SEqHatDot_4);
    SEqHatDot_1 = SEqHatDot_1 / norm;
    SEqHatDot_2 = SEqHatDot_2 / norm;
    SEqHatDot_3 = SEqHatDot_3 / norm;
    SEqHatDot_4 = SEqHatDot_4 / norm;

// compute angular estimated direction of the gyroscope error
    w_err_x = twoSEq_1 * SEqHatDot_2 - twoSEq_2 * SEqHatDot_1 - twoSEq_3 * SEqHatDot_4 + twoSEq_4 * SEqHatDot_3;
    w_err_y = twoSEq_1 * SEqHatDot_3 + twoSEq_2 * SEqHatDot_4 - twoSEq_3 * SEqHatDot_1 - twoSEq_4 * SEqHatDot_2;
    w_err_z = twoSEq_1 * SEqHatDot_4 - twoSEq_2 * SEqHatDot_3 + twoSEq_3 * SEqHatDot_2 - twoSEq_4 * SEqHatDot_1;

// compute and remove the gyroscope baises
    w_bx += w_err_x * deltat * zeta;
    w_by += w_err_y * deltat * zeta;
    w_bz += w_err_z * deltat * zeta;

    w_x -= w_bx;
    w_y -= w_by;
    w_z -= w_bz;
    //陀螺仪偏移减掉

// compute the quaternion rate measured by gyroscopes
    SEqDot_omega_1 = -halfSEq_2 * w_x - halfSEq_3 * w_y - halfSEq_4 * w_z;
    SEqDot_omega_2 =  halfSEq_1 * w_x + halfSEq_3 * w_z - halfSEq_4 * w_y;
    SEqDot_omega_3 =  halfSEq_1 * w_y - halfSEq_2 * w_z + halfSEq_4 * w_x;
    SEqDot_omega_4 =  halfSEq_1 * w_z + halfSEq_2 * w_y - halfSEq_3 * w_x;

// compute then integrate the estimated quaternion rate
    SEq_1 += (SEqDot_omega_1 - (beta * SEqHatDot_1)) * deltat;
    SEq_2 += (SEqDot_omega_2 - (beta * SEqHatDot_2)) * deltat;
    SEq_3 += (SEqDot_omega_3 - (beta * SEqHatDot_3)) * deltat;
    SEq_4 += (SEqDot_omega_4 - (beta * SEqHatDot_4)) * deltat;
    //姿态融合，梯度下降

// normalise quaternion
    norm = sqrt(SEq_1 * SEq_1 + SEq_2 * SEq_2 + SEq_3 * SEq_3 + SEq_4 * SEq_4);
    SEq_1 /= norm;
    SEq_2 /= norm;
    SEq_3 /= norm;
    SEq_4 /= norm;

// compute flux in the earth frame
    SEq_1SEq_2 = SEq_1 * SEq_2;

// recompute axulirary variables
    SEq_1SEq_3 = SEq_1 * SEq_3;
    SEq_1SEq_4 = SEq_1 * SEq_4;
    SEq_3SEq_4 = SEq_3 * SEq_4;
    SEq_2SEq_3 = SEq_2 * SEq_3;
    SEq_2SEq_4 = SEq_2 * SEq_4;

    h_x = twom_x * (0.5f - SEq_3 * SEq_3 - SEq_4 * SEq_4) + twom_y * (SEq_2SEq_3 - SEq_1SEq_4) + twom_z * (SEq_2SEq_4 + SEq_1SEq_3);
    h_y = twom_x * (SEq_2SEq_3 + SEq_1SEq_4) + twom_y * (0.5f - SEq_2 * SEq_2 - SEq_4 * SEq_4) + twom_z * (SEq_3SEq_4 - SEq_1SEq_2);
    h_z = twom_x * (SEq_2SEq_4 - SEq_1SEq_3) + twom_y * (SEq_3SEq_4 + SEq_1SEq_2) + twom_z * (0.5f - SEq_2 * SEq_2 - SEq_3 * SEq_3);

// normalise the flux vector to have only components in the x and z
    b_x = sqrt((h_x * h_x) + (h_y * h_y));
    b_z = h_z;
}


SPISettings bmi323Settings(SPI_FREQ, MSBFIRST, SPI_MODE0);
unsigned long lastPrintTime = 0;
const unsigned long printInterval = 1000 / OUTPUT_RATE_HZ;

uint16_t readRegister16(SPIClass &spi, uint8_t cs, uint8_t reg) {
  spi.beginTransaction(bmi323Settings);
  digitalWrite(cs, LOW);
  delayMicroseconds(1);

  spi.transfer(reg | 0x80);
  spi.transfer(0x00);
  uint8_t lsb = spi.transfer(0x00);
  uint8_t msb = spi.transfer(0x00);

  digitalWrite(cs, HIGH);
  spi.endTransaction();
  delayMicroseconds(2);

  return (msb << 8) | lsb;
}

void writeRegister16(SPIClass &spi, uint8_t cs, uint8_t reg, uint16_t data) {
  spi.beginTransaction(bmi323Settings);
  digitalWrite(cs, LOW);
  delayMicroseconds(1);

  spi.transfer(reg & 0x7F);
  spi.transfer(data & 0xFF);
  spi.transfer((data >> 8) & 0xFF);

  digitalWrite(cs, HIGH);
  spi.endTransaction();
  delayMicroseconds(2);
}

bool initializeFeatureEngine(SPIClass &spi, uint8_t cs) {
  writeRegister16(spi, cs, FEATURE_IO0_REG, 0x0000);
  delay(1);

  writeRegister16(spi, cs, FEATURE_IO2_REG, 0x012C);
  delay(1);

  writeRegister16(spi, cs, FEATURE_IO_STATUS_REG, 0x0001);
  delay(1);

  writeRegister16(spi, cs, FEATURE_CTRL_REG, 0x0001);
  delay(10);

  int timeout = 0;
  uint16_t featureIO1Status;
  do {
    delay(10);
    featureIO1Status = readRegister16(spi, cs, FEATURE_IO1_REG);
    uint8_t errorStatus = featureIO1Status & 0x0F;
    if (errorStatus == 0x01) return true;
    if (errorStatus == 0x03) return false;
    timeout++;
  } while ((featureIO1Status & 0x0F) == 0x00 && timeout < 50);

  return (timeout < 50);
}

float convertAccelData(uint16_t rawData) {
  int16_t signedData = (int16_t)rawData;
  return signedData / ACC_RANGE_LSB_PER_G;
}
float convertGyroData(uint16_t rawData) {
  int16_t signedData = (int16_t)rawData;
  return signedData / GYR_RANGE_LSB_PER_DPS;
}
float convertTempData(uint16_t rawData) {
  int16_t signedData = (int16_t)rawData;
  return (signedData / 512.0) + 23.0;
}

bool initBMI323(SPIClass &spi, uint8_t cs) {
  readRegister16(spi, cs, CHIP_ID_REG);
  delayMicroseconds(200);

  uint16_t chipID = readRegister16(spi, cs, CHIP_ID_REG);
  Serial.print("Chip ID: 0x");
  Serial.println(chipID, HEX);

  if ((chipID & 0xFF) != 0x43 && (chipID & 0xFF) != 0x41) {
    return false;
  }

  writeRegister16(spi, cs, CMD_REG, SOFT_RESET_CMD);
  delay(5);

  if (!initializeFeatureEngine(spi, cs)) return false;

  writeRegister16(spi, cs, ACC_CONF_REG, ACC_CONF_NORMAL_HZ_G);
  writeRegister16(spi, cs, GYR_CONF_REG, GYR_CONF_NORMAL_HZ_DPS);
  delay(50);

  return true;
}

void setup() {
  Serial.begin(115200);
  delay(2000);

  //BMI1
  pinMode(BMI1_CS, OUTPUT);
  digitalWrite(BMI1_CS, HIGH);
  spi1.begin(BMI1_SCK, BMI1_MISO, BMI1_MOSI, BMI1_CS);

  //BMI2
  pinMode(BMI2_CS, OUTPUT);
  digitalWrite(BMI2_CS, HIGH);
  spi2.begin(BMI2_SCK, BMI2_MISO, BMI2_MOSI, BMI2_CS);

  delay(100);

  Serial.println("Initializing BMI1...");
  bool ok1 = initBMI323(spi1, BMI1_CS);
  Serial.println(ok1 ? "BMI1 Ready" : "BMI1 Init Failed!");

  Serial.println("Initializing BMI2...");
  bool ok2 = initBMI323(spi2, BMI2_CS);
  Serial.println(ok2 ? "BMI2 Ready" : "BMI2 Init Failed!");

  if (!ok1 && !ok2) {
    Serial.println("No BMI323 found, check wiring!");
    while (1);
  }

  Serial.println("Format: ID, ax,ay,az,gx,gy,gz,temp");
}

void loop() {
  unsigned long currentTime = millis();
  if (currentTime - lastPrintTime < printInterval) return;
  lastPrintTime = currentTime;
  //BMI1
  uint16_t accX1 = readRegister16(spi1, BMI1_CS, ACC_DATA_X_REG);
  uint16_t accY1 = readRegister16(spi1, BMI1_CS, ACC_DATA_Y_REG);
  uint16_t accZ1 = readRegister16(spi1, BMI1_CS, ACC_DATA_Z_REG);
  uint16_t gyrX1 = readRegister16(spi1, BMI1_CS, GYR_DATA_X_REG);
  uint16_t gyrY1 = readRegister16(spi1, BMI1_CS, GYR_DATA_Y_REG);
  uint16_t gyrZ1 = readRegister16(spi1, BMI1_CS, GYR_DATA_Z_REG);
  uint16_t tempRaw1 = readRegister16(spi1, BMI1_CS, TEMP_DATA_REG);

  //BMI2
  uint16_t accX2 = readRegister16(spi2, BMI2_CS, ACC_DATA_X_REG);
  uint16_t accY2 = readRegister16(spi2, BMI2_CS, ACC_DATA_Y_REG);
  uint16_t accZ2 = readRegister16(spi2, BMI2_CS, ACC_DATA_Z_REG);
  uint16_t gyrX2 = readRegister16(spi2, BMI2_CS, GYR_DATA_X_REG);
  uint16_t gyrY2 = readRegister16(spi2, BMI2_CS, GYR_DATA_Y_REG);
  uint16_t gyrZ2 = readRegister16(spi2, BMI2_CS, GYR_DATA_Z_REG);
  uint16_t tempRaw2 = readRegister16(spi2, BMI2_CS, TEMP_DATA_REG);
/*
  //BMI1
  Serial.print("BMI1,");
  Serial.print(convertAccelData(accX1), 3); Serial.print(",");
  Serial.print(convertAccelData(accY1), 3); Serial.print(",");
  Serial.print(convertAccelData(accZ1), 3); Serial.print(",");
  Serial.print(convertGyroData(gyrX1), 2); Serial.print(",");
  Serial.print(convertGyroData(gyrY1), 2); Serial.print(",");
  Serial.print(convertGyroData(gyrZ1), 2); Serial.print(",");
  Serial.println(convertTempData(tempRaw1), 1);
  
  //BMI2
  Serial.print("BMI2,");
  Serial.print(convertAccelData(accX2), 3); Serial.print(",");
  Serial.print(convertAccelData(accY2), 3); Serial.print(",");
  Serial.print(convertAccelData(accZ2), 3); Serial.print(",");
  Serial.print(convertGyroData(gyrX2), 2); Serial.print(",");
  Serial.print(convertGyroData(gyrY2), 2); Serial.print(",");
  Serial.print(convertGyroData(gyrZ2), 2); Serial.print(",");
  Serial.println(convertTempData(tempRaw2), 1);*/
  filterUpdate(convertGyroData(gyrX1)*Gyro_Gr, convertGyroData(gyrY1)*Gyro_Gr, convertGyroData(gyrZ1)*Gyro_Gr, convertAccelData(accX1), convertAccelData(accY1), convertAccelData(accZ1), 0, 0, 0);
  Serial.print("BMI1:");Serial.print(SEq_1);Serial.print(",");Serial.print(SEq_2);Serial.print(",");Serial.print(SEq_3);Serial.print(",");Serial.println(SEq_4);
}
