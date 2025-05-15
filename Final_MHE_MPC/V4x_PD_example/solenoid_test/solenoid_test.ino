/*******************************************************************************
 * Maggy V4.x - Solenoid Sequential Test (Modified)
 * 
 * Modified Test Sequence:
 * 1. Initialize the system normally (sensors, solenoids, calibration)
 * 2. Apply a control signal of 100 to X+ solenoid for 1 period
 * 3. No control signal for 1 period
 * 4. Apply a control signal of -100 to X+ solenoid for 1 period
 * 5. No control signal for 1 period
 * 6. Apply a control signal of 100 to X- solenoid for 1 period
 * 7. No control signal for 1 period
 * 8. Apply a control signal of -100 to X- solenoid for 1 period
 * 9. No control signal for 1 period
 * 10. Apply a control signal of 100 to Y+ solenoid for 1 period
 * 11. No control signal for 1 period
 * 12. Apply a control signal of -100 to Y+ solenoid for 1 period
 * 13. No control signal for 1 period
 * 14. Apply a control signal of 100 to Y- solenoid for 1 period
 * 15. No control signal for 1 period
 * 16. Apply a control signal of -100 to Y- solenoid for 1 period
 * 17. No control signal for 1 period
 * 18. Apply a control signal of 100 to all solenoids for 1 period
 * 19. No control signal for 1 period
 * 20. Apply a control signal of -100 to all solenoids for 1 period
 * 21. No control signal for 1 period
 * 22. Apply a control signal of 100 to X+ and -100 to X- for 1 period
 * 23. No control signal for 1 period
 * 24. Apply a control signal of 100 to Y+ and -100 to Y- for 1 period
 * 25. No control signal for 3 periods
 *
 * The period duration is adjustable via the PERIOD_DURATION constant.
 *******************************************************************************/

#include "definitions.h"
#include "functions.h"

// Sensor configuration - channels defined here; NUM_SENSORS and PRIMARY_SENSOR_INDEX in definitions.h
const int SENSOR_CHANNELS[NUM_SENSORS] = {7, 2, 3}; // Sensor channel IDs on the multiplexer (also indicated on the PCB)

// Test parameters
constexpr int TEST_SIGNAL_AMPLITUDE = 100;    // Fixed control signal amplitude
constexpr unsigned long PERIOD_DURATION = 1000000*0.2;  // 1 period = 1,000,000 microseconds (i.e. 1 second)
// (Change PERIOD_DURATION to adjust the period length)

// Structure defining a test step
struct TestStep {
  float xPos;
  float xNeg;
  float yPos;
  float yNeg;
  unsigned long duration;  // Duration of this step in microseconds
};

// Test sequence array definition (each step uses periods as defined above)
constexpr int NUM_TEST_STEPS = 24; // There are 24 transitions; note that the last step uses 3 periods.
TestStep testSequence[NUM_TEST_STEPS + 1] = {
  // Step 2: 100 to X+ for 1 period
  { TEST_SIGNAL_AMPLITUDE, 0, 0, 0, PERIOD_DURATION },
  // Step 3: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 4: -100 to X+ for 1 period
  { -TEST_SIGNAL_AMPLITUDE, 0, 0, 0, PERIOD_DURATION },
  // Step 5: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 6: 100 to X- for 1 period
  { 0, TEST_SIGNAL_AMPLITUDE, 0, 0, PERIOD_DURATION },
  // Step 7: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 8: -100 to X- for 1 period
  { 0, -TEST_SIGNAL_AMPLITUDE, 0, 0, PERIOD_DURATION },
  // Step 9: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 10: 100 to Y+ for 1 period
  { 0, 0, TEST_SIGNAL_AMPLITUDE, 0, PERIOD_DURATION },
  // Step 11: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 12: -100 to Y+ for 1 period
  { 0, 0, -TEST_SIGNAL_AMPLITUDE, 0, PERIOD_DURATION },
  // Step 13: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 14: 100 to Y- for 1 period
  { 0, 0, 0, TEST_SIGNAL_AMPLITUDE, PERIOD_DURATION },
  // Step 15: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 16: -100 to Y- for 1 period
  { 0, 0, 0, -TEST_SIGNAL_AMPLITUDE, PERIOD_DURATION },
  // Step 17: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 18: 100 to all solenoids for 1 period
  { TEST_SIGNAL_AMPLITUDE, TEST_SIGNAL_AMPLITUDE, TEST_SIGNAL_AMPLITUDE, TEST_SIGNAL_AMPLITUDE, PERIOD_DURATION },
  // Step 19: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 20: -100 to all solenoids for 1 period
  { -TEST_SIGNAL_AMPLITUDE, -TEST_SIGNAL_AMPLITUDE, -TEST_SIGNAL_AMPLITUDE, -TEST_SIGNAL_AMPLITUDE, PERIOD_DURATION },
  // Step 21: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 22: 100 to X+ and -100 to X- for 1 period
  { TEST_SIGNAL_AMPLITUDE, -TEST_SIGNAL_AMPLITUDE, 0, 0, PERIOD_DURATION },
  // Step 23: No signal for 1 period
  { 0, 0, 0, 0, PERIOD_DURATION },
  // Step 24: 100 to Y+ and -100 to Y- for 1 period
  { 0, 0, TEST_SIGNAL_AMPLITUDE, -TEST_SIGNAL_AMPLITUDE, PERIOD_DURATION },
  // Step 25: No signal for 3 periods
  { 0, 0, 0, 0, 3*PERIOD_DURATION }
};

// Sensor objects - one per sensor
Tlv493d Sensors[NUM_SENSORS];
TCA9548 mux_sensors(0x70);

// Timing parameters for sensor reading and logging
constexpr float sensorFrequency = 5000.0;
constexpr int sensorInterval = round(1e6 / sensorFrequency);
constexpr float loggingFrequency = 500.0;    // Log at 100 Hz
constexpr int loggingInterval = round(1e6 / loggingFrequency);

// Timing variables
unsigned long prevSensorTime = 0;
unsigned long prevLoggingTime = 0;
unsigned long testStepStartTime = 0;
float realSamplingFreq = 0;
int loggingCounter = 0;

// Control variables for solenoid currents (and sensor values)
float magFieldX = 0, magFieldY = 0, magFieldZ = 0;
float prevMagFieldX = 0, prevMagFieldY = 0, prevMagFieldZ = 0;
float currentXPos = 0, currentXNeg = 0, currentYPos = 0, currentYNeg = 0;
float pwmInputX = 0, pwmInputY = 0;

// Variables for sensor freeze detection
float lastMagField[NUM_SENSORS][3] = {{0}};
int freezeCounter = 0;

// Global index to keep track of the current test sequence step
int testStepIndex = 0;

void setup(){
  Serial.begin(115200);
  delay(100);
  Serial.println("Solenoid Sequential Test - Starting initialization");

  initializeSensors();
  delay(50);

  initializeSolenoids();
  delay(50);

  calibrateSensors();
  delay(50);

  // calibrateDirectFeedthrough(); // Uncoment to test if feedthrough is removed correctly
  delay(50);

  // Other initializations if needed...
  delay(50);
  
  Serial.println("Initialization complete");
  Serial.println("Starting modified solenoid test sequence");

  // Start test sequence:
  testStepIndex = 0;
  testStepStartTime = micros();
  // Ensure all solenoids start at 0:
  applyTestSignals(0, 0, 0, 0);
  delay(1000); // Brief pause before starting test

  // Apply signals for the first test step:
  updateControlSignals();
}

void loop(){
  unsigned long currentTime = micros();

  // Sensor reading and processing (same as in main program)
  if(currentTime - prevSensorTime >= (unsigned long)sensorInterval){
    unsigned long dt = currentTime - prevSensorTime;
    if(dt > 0) realSamplingFreq = 1e6 / (float)dt;

    currentXPos = getSolenoidCurrent(CURRENT_X_POS);
    currentXNeg = getSolenoidCurrent(CURRENT_X_NEG);
    currentYPos = getSolenoidCurrent(CURRENT_Y_POS);
    currentYNeg = getSolenoidCurrent(CURRENT_Y_NEG);

    readAllSensors();
    processSensorData(currentXPos, currentXNeg, currentYPos, currentYNeg);

    magFieldX = rawMagFieldDetrended[PRIMARY_SENSOR_INDEX][0];
    magFieldY = rawMagFieldDetrended[PRIMARY_SENSOR_INDEX][1];
    magFieldZ = rawMagFieldDetrended[PRIMARY_SENSOR_INDEX][2];

    prevMagFieldX = magFieldX;
    prevMagFieldY = magFieldY;
    prevMagFieldZ = magFieldZ;
    prevSensorTime = currentTime;

    checkSensorFreeze(lastMagField, &freezeCounter);
  }

  // Test sequence state machine:
  unsigned long timeInStep = currentTime - testStepStartTime;
  if (timeInStep >= testSequence[testStepIndex].duration) {
    // Move to the next step
    testStepIndex++;
    if (testStepIndex >= (NUM_TEST_STEPS + 1)) { 
      // Option: Restart the test sequence continuously.
      testStepIndex = 0;
    }
    testStepStartTime = currentTime; // Restart timer for the new step

    updateControlSignals();
  }

  // Data logging (occurs at the specified logging frequency)
  if (currentTime - prevLoggingTime >= (unsigned long)loggingInterval) {
    logSystemState(currentTime, currentXPos, currentXNeg, currentYPos, currentYNeg);
    loggingCounter++;
    prevLoggingTime = currentTime;
  }
}

// Update and apply solenoid control signals for the current test step
void updateControlSignals() {
  // Get current step values from the test sequence array
  float xPos = testSequence[testStepIndex].xPos;
  float xNeg = testSequence[testStepIndex].xNeg;   
  float yPos = testSequence[testStepIndex].yPos;
  float yNeg = testSequence[testStepIndex].yNeg;

  applyTestSignals(xPos, xNeg, yPos, yNeg);
}

// Apply individual control signals to each solenoid
void applyTestSignals(float xPos, float xNeg, float yPos, float yNeg) {
  setSolenoidInput(xPos, MD2_IN1, MD2_IN2);   // X+
  setSolenoidInput(xNeg, MD3_IN1, MD3_IN2);   // X-
  setSolenoidInput(yPos, MD4_IN1, MD4_IN2);   // Y+
  setSolenoidInput(yNeg, MD1_IN1, MD1_IN2);   // Y-
}
