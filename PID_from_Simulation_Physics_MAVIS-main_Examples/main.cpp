#include "../../src/interface/interface.h"

struct PID_Controller
{
    double kp;
    double ki;
    double kd;
    double previous_error; // see if setting to 0.0 is necessary
    double integral;
    double controllerP;
    double controllerI;
    double controllerD;
    double clegg; // clegg controller for integral term (not used right now)
    /**
     * get_errors Calculates the errors and adjusts the system using the PID controller algorithm.
     *
     * @param kp The proportional gain.
     * @param ki The integral gain.
     * @param kd The derivative gain.
     * @param setpoint The desired setpoint value.
     * @param actual_value The actual value of the system.
     * @param dt The time step.
     */
    void get_errors(double kp, double ki, double kd, double setpoint, double actual_value, double dt) {
        double error = setpoint - actual_value;
        controllerP = kp * error;
        // If error has changed sign, reset integral to prevent windup
        // if (previous_error * error < 0) {
        //     integral = 0;
        // } else {
        //     integral += 0.5 * (error + previous_error) * dt;  // Trapezoidal rule
        // }
        integral += 0.5 * (error + previous_error) * dt;
        if (integral > 0.5)  // Anti-windup for integral term
            integral = 0.5;
        else if (integral < -0.5)
            integral = -0.5;
        controllerI = ki * integral;
        controllerD = kd * (error - previous_error) / dt;
        previous_error = error;  // Update previous error
    }

    /**
     * get_command Calculates the output command based on the pilot control, commanded rate, and controller parameters.
     * 
     * @param pilotControl The pilot control input.
     * @param commandedRate The commanded rate input.
     * @param output The calculated output command.
     */
    void get_command(double pilotControl, double commandedRate, bool gettingFinalCommand, double& output) {
        if (gettingFinalCommand == true) // fix this logic!!!
            {
            if (pilotControl <= 0.01 && pilotControl >= -0.01)
            {
                output = controllerP + controllerI + controllerD;
            }
            else
            {
                output = pilotControl;
            }
            }
        else
            {
                if (pilotControl <= 0.01 && pilotControl >= -0.01)
                    {
                        output = controllerP + controllerI + controllerD;
                    }
                    else
                    {
                    output = commandedRate;
                    }
            }
    }
}; // End of PID_Controller struct (you have to end structs with a semicolon)

int main (int argc, char * const argv[]) {
    if (argc != 2) {
        cout<<"Please enter a filename on the command line in the form: " << argv[0] << " <filename.json>" << endl;
        return 1;
    }
    // Initialize variables
    double d_error_p = 0.0;
    double old_time = 0.0;
    double dt = 0.0;
    double p_command = 0.0;
    double q_command = 0.0;
    double r_command = 0.0;
    double phi_command = 0.0;
    double theta_command = 0.0;
    double psi_command = 0.0;
    double u_command = 0.0;
    double alpha = 0.0;

    const char* filename = argv[1];
    cout << "Program started with file: " << filename << endl;

    ifstream f(filename);
    json inputDict = json::parse(f);

    double endTime = json_get_double(inputDict,"end_time[s]");
    json interfaceDict = json_get_dictionary(inputDict,"interface");

    Interface* myInterface = new Interface();

    // Add connections(s)
    string name = "physics_to_controller";
    json connectionDict = json_get_dictionary(interfaceDict,name);
    Connection* mPhysicsToController = new Connection(connectionDict,name);
    myInterface->add_connection(mPhysicsToController);

    name = "controller_to_physics";
    connectionDict = json_get_dictionary(interfaceDict,name);
    Connection* mControllerToPhysics = new Connection(connectionDict,name);
    myInterface->add_connection(mControllerToPhysics);

    name = "pilot_to_controller";
    connectionDict = json_get_dictionary(interfaceDict,name);
    Connection* mPilotToController = new Connection(connectionDict,name);
    myInterface->add_connection(mPilotToController);

    // Set up PID controllers
    PID_Controller aileronPID; // directly control roll rate using aileron
    json aileronDict = json_get_dictionary(inputDict,"aileron");
    p_command = json_get_double(aileronDict,"rollRateCommand");
    aileronPID.kp = json_get_double(aileronDict,"P");
    aileronPID.ki = json_get_double(aileronDict,"I");
    aileronPID.kd = json_get_double(aileronDict,"D");

    PID_Controller phiPID; // control roll angle using aileron (cascade with roll rate controller)
    json phiDict = json_get_dictionary(inputDict,"phi");
    phi_command = json_get_double(phiDict,"rollCommand");
    phiPID.kp = json_get_double(phiDict,"P");
    phiPID.ki = json_get_double(phiDict,"I");
    phiPID.kd = json_get_double(phiDict,"D");

    PID_Controller elevatorPID; // control pitch rate using elevator
    json elevatorDict = json_get_dictionary(inputDict,"elevator");
    q_command = json_get_double(elevatorDict,"pitchRateCommand");
    elevatorPID.kp = json_get_double(elevatorDict,"P");
    elevatorPID.ki = json_get_double(elevatorDict,"I");
    elevatorPID.kd = json_get_double(elevatorDict,"D");

    PID_Controller thetaPID; // control pitch angle using elevator (cascade with pitch rate controller)
    json thetaDict = json_get_dictionary(inputDict,"theta");
    theta_command = json_get_double(thetaDict,"pitchCommand");
    thetaPID.kp = json_get_double(thetaDict,"P");
    thetaPID.ki = json_get_double(thetaDict,"I");
    thetaPID.kd = json_get_double(thetaDict,"D");

    PID_Controller rudderPID; // control yaw rate using rudder
    json rudderDict = json_get_dictionary(inputDict,"rudder");
    r_command = json_get_double(rudderDict,"yawRateCommand");
    rudderPID.kp = json_get_double(rudderDict,"P");
    rudderPID.ki = json_get_double(rudderDict,"I");
    rudderPID.kd = json_get_double(rudderDict,"D");

    PID_Controller throttlePID; // control velocity using throttle
    json throttleDict = json_get_dictionary(inputDict,"throttle");
    u_command = json_get_double(throttleDict,"u_command");
    throttlePID.kp = json_get_double(throttleDict,"P");
    throttlePID.ki = json_get_double(throttleDict,"I");
    throttlePID.kd = json_get_double(throttleDict,"D");

    // Run until initialized states from the physics are received
    printf("Waiting for initial data from sim...\n");
    do {
        myInterface->read(); // only read until we establish a connection with the sim
    } while(mPhysicsToController->mArrayValues[0] == 0.0); // note it was initialized to zero automatically.

    // Copy current states of the aircraft into local working variables
    double simTime = mPhysicsToController->mArrayValues[0];

    double u = mPhysicsToController->mArrayValues[1];
    double v = mPhysicsToController->mArrayValues[2];
    double w = mPhysicsToController->mArrayValues[3];

    double p = mPhysicsToController->mArrayValues[4];
    double q = mPhysicsToController->mArrayValues[5];
    double r = mPhysicsToController->mArrayValues[6];

    double xf = mPhysicsToController->mArrayValues[7];
    double yf = mPhysicsToController->mArrayValues[8];
    double zf = mPhysicsToController->mArrayValues[9];

    double phi   = mPhysicsToController->mArrayValues[10];
    double theta = mPhysicsToController->mArrayValues[11];
    double psi   = mPhysicsToController->mArrayValues[12];

    double da  = mPhysicsToController->mArrayValues[13];
    double de  = mPhysicsToController->mArrayValues[14];
    double dr  = mPhysicsToController->mArrayValues[15];
    double tau = mPhysicsToController->mArrayValues[16];

    double dapilot  = mPilotToController->mArrayValues[0];
    double depilot  = mPilotToController->mArrayValues[1];
    double drpilot  = mPilotToController->mArrayValues[2];
    double taupilot = mPilotToController->mArrayValues[3];

    printf("Initial data read from Sim.\n");
    printf("Initial data read from Sim.\n");
    printf("Sim time = %f\n",simTime);
    printf("Sim throttle = %f\n",tau);
    
    int loopCount = 0; // Counter for the number of loops
    do {    
            // Send/Receive
            myInterface->update(); // Updates the interface. This only sends/receives if it is time to send based on the update hz specified in the json file

            // If there is new data to be used
            if(mPhysicsToController->mArrayValues[0] > simTime) 
            {
                old_time = simTime; // Stores the old simulation time
                // Copy out new states of the aircraft as needed for controller
                simTime = mPhysicsToController->mArrayValues[0]; // Current simulation time
                p = mPhysicsToController->mArrayValues[4]; // Roll rate
                q = mPhysicsToController->mArrayValues[5]; // Pitch rate
                r = mPhysicsToController->mArrayValues[6]; // Yaw rate
                phi = mPhysicsToController->mArrayValues[10]; // Roll angle
                theta = mPhysicsToController->mArrayValues[11]; // Pitch angle
                tau = mPhysicsToController->mArrayValues[16]; // Throttle command
                u = mPhysicsToController->mArrayValues[1]; // Body-axis x velocity
                v = mPhysicsToController->mArrayValues[2]; // Body-axis y velocity
                w = mPhysicsToController->mArrayValues[3]; // Body-axis z velocity
                dt = simTime - old_time; // Time difference

                dapilot = pow(mPilotToController->mArrayValues[0], 3); // Aileron control input (using exponents so that it's not too sensitive)
                depilot = pow(mPilotToController->mArrayValues[1], 3); // Elevator control input
                drpilot = pow(mPilotToController->mArrayValues[2], 3); // Rudder control input
                taupilot = mPilotToController->mArrayValues[3]; // Throttle control input

                alpha = atan2(w, u); // Angle of attack
                theta_command = alpha; // Set point for pitch angle

                // cascade portion. Calculate errors for roll and pitch
                phiPID.get_errors(phiPID.kp, phiPID.ki, phiPID.kd, phi_command, phi, dt); // Calculate PID errors for roll angle (phi)
                thetaPID.get_errors(thetaPID.kp, thetaPID.ki, thetaPID.kd, theta_command, theta, dt); // Calculate PID errors for pitch angle (theta)
                phiPID.get_command(dapilot, p_command, false, p_command); // getting roll rate (p) command after controller calculations
                thetaPID.get_command(depilot, q_command, false, q_command); // getting pitch rate (q) command after controller calculations

                // Calculate errors for p, q, r, and u, notice that the setpoints are the commands from the cascaded controllers
                aileronPID.get_errors(aileronPID.kp, aileronPID.ki, aileronPID.kd, p_command, p, dt); // Calculate PID errors for roll rate (p)
                elevatorPID.get_errors(elevatorPID.kp, elevatorPID.ki, elevatorPID.kd, q_command, q, dt); // Calculate PID errors for pitch rate (q)
                rudderPID.get_errors(rudderPID.kp, rudderPID.ki, rudderPID.kd, r_command, r, dt); // Calculate PID errors for yaw rate (r)
                throttlePID.get_errors(throttlePID.kp, throttlePID.ki, throttlePID.kd, u_command, u, dt); // Calculate PID errors for velocity (u)
                
                // Calculate final commands based on the cascaded and inner controllers
                throttlePID.get_command(taupilot, tau, true, tau); // getting final throttle command after controller calculations
                aileronPID.get_command(dapilot, da, true, da); // getting final aileron command after controller calculations      
                elevatorPID.get_command(depilot, de, true, de); // getting final elevator command after controller calculations
                rudderPID.get_command(drpilot, dr, true, dr); // getting final rudder command after controller calculations

                // Set commanded controls to send to the physics
                mControllerToPhysics->mArrayValues[0] = da; // Aileron control command
                mControllerToPhysics->mArrayValues[1] = de; // Elevator control command
                mControllerToPhysics->mArrayValues[2] = dr; // Rudder control command
                mControllerToPhysics->mArrayValues[3] = tau; // Throttle control command
                loopCount++; // Increment loop counter
            
                printf("Sim time = %f, aileron: %f\n", simTime, da);
                printf("Sim time = %f, elevator: %f\n", simTime, de);
                printf("Sim time = %f, rudder: %f\n", simTime, dr);
                printf("Sim time = %f, p: %f\n", simTime, p);
                printf("Sim time = %f, q: %f\n", simTime, q);
                printf("Sim time = %f, r: %f\n", simTime, r);
                printf("Sim time = %f, roll angle: %f\n", simTime, phi);
                printf("Sim time = %f, pitch angle: %f\n", simTime, theta);
                printf("Sim time = %f, alpha: %f\n", simTime, alpha);
                printf("Sim time = %f, yaw angle: %f\n", simTime, psi);
                printf("Sim time = %f, u: %f\n", simTime, u);
                printf("Sim time = %f, throttle: %f\n", simTime, tau);
                // printf("Sim time = %f, Aileron PID is rate controller: %d\n", simTime, aileronPID.isRateController);
                // printf("Sim time = %f, Elevator PID is rate controller: %d\n", simTime, elevatorPID.isRateController);
                // printf("Sim time = %f, Rudder PID is rate controller: %d\n", simTime, rudderPID.isRateController);
                // printf("Sim time = %f, Throttle PID is rate controller: %d\n", simTime, throttlePID.isRateController);
                // printf("Sim time = %f, Theta PID is rate controller: %d\n", simTime, thetaPID.isRateController);
            }
    } while (simTime < endTime); // Continue until the simulation time is less than the end time
}
