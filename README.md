# Orientation Filter

This project implements an orientation filter for a MARG (Magnetometer, Accelerometer, and Gyroscope) sensor array in Swift. The filter compensates for magnetic distortion and gyroscope drift, ensuring accurate orientation estimation.

## Features

- **Gyroscope Drift Compensation**: Reduces errors caused by gyroscope drift.
- **Magnetic Distortion Compensation**: Adjusts for magnetic interference to provide accurate orientation data.
- **Quaternion-Based Orientation Estimation**: Uses quaternions for smooth and continuous orientation representation.

## Requirements

- Swift 2.0 or later

## Installation

Include the `Filter` class in your Swift project.

## Usage

To use the orientation filter, create an instance of the `Filter` class and update it with sensor readings.

### Example

```swift
import Foundation

// Create an instance of the Filter class
let filter = Filter()

// Update the filter with current sensor readings
filter.filterUpdate(
    w_x: gyroX, w_y: gyroY, w_z: gyroZ, 
    a_x: accelX, a_y: accelY, a_z: accelZ, 
    m_x: magX, m_y: magY, m_z: magZ
)
```

### Parameters

- `w_x`, `w_y`, `w_z`: Gyroscope readings (in radians per second).
- `a_x`, `a_y`, `a_z`: Accelerometer readings.
- `m_x`, `m_y`, `m_z`: Magnetometer readings.

### How It Works

1. **Normalization**: The accelerometer and magnetometer readings are normalized.
2. **Objective Function and Jacobian**: Computes the objective function and its Jacobian to estimate the direction of the gyroscope error.
3. **Gradient Descent**: The filter uses gradient descent to minimize the error in the orientation estimate.
4. **Quaternion Update**: Updates the quaternion representing the current orientation.
5. **Gyroscope Bias Correction**: Corrects for any bias in the gyroscope readings.

## Contributing

Contributions are welcome. Please submit a pull request or open an issue to discuss any changes.

## Author

- Adam Leitgeb  
- Petr Chmelař

## License

This project is licensed under the MIT License. See the LICENSE file for details.
