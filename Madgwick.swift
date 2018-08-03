import Foundation

public final class Filter {

    // MARK: System constants

    /// sampling period in seconds (shown as 1ms)
    private let deltat = 0.001

    /// gyroscope measurement error in rad/s (shown as 5 deg/s)
    private static let gyroMassError = 3.14159265358979 * (5 / 180)

    /// gyroscope measurement error in rad/s/s (shown as 0.2 deg/s/s)
    private static let gyroMassDrift = 3.14159265358979 * (2 / 180)

    // compute beta
    private let beta = sqrt(3.0 / 4.0) * Madgwick.gyroMassError

    // compute zeta
    private let zeta = sqrt(3.0 / 4.0) * Madgwick.gyroMassDrift

    // MARK: Global system variables

    // estimated orientation quaternion elements with initial conditions
    var SEq_1: Double = 1.0
    var SEq_2: Double = 0.0
    var SEq_3: Double = 0.0
    var SEq_4: Double = 0.0

    // reference direction o flux in earth frame
    var b_x: Double = 1.0
    var b_z: Double = 0.0

    // estimate gyroscope biases error
    var w_bx: Double = 0.0
    var w_by: Double = 0.0
    var w_bz: Double = 0.0

    // MARK: - Methods

    /// method to compute one filter iteration
    ///
    /// - Parameters:
    ///   - w_x: gyroscope x value in rad/s
    ///   - w_y: gyroscope y value in rad/s
    ///   - w_z: gyroscope z value in rad/s
    ///   - a_x: accelerometer x value
    ///   - a_y: accelerometer y value
    ///   - a_z: accelerometer z value
    ///   - m_x: magnetometer x value
    ///   - m_y: magnetometer y value
    ///   - m_z: magnetometer z value
    public func filterUpdate(w_x: Double, w_y: Double, w_z: Double, a_x: Double, a_y: Double, a_z: Double, m_x: Double, m_y: Double, m_z: Double) {

        var mutable_w_x: Double = w_x
        var mutable_w_y: Double = w_y
        var mutable_w_z: Double = w_z

        var mutable_a_x: Double = a_x
        var mutable_a_y: Double = a_z
        var mutable_a_z: Double = a_z

        var mutable_m_x: Double = m_x
        var mutable_m_y: Double = m_y
        var mutable_m_z: Double = m_z

        // vector norm

        var norm: Double!

        // quaternion rate from gyroscopes elements

        var SEqDot_omega_1: Double!
        var SEqDot_omega_2: Double!
        var SEqDot_omega_3: Double!
        var SEqDot_omega_4: Double!

        // objective function elements

        var f_1: Double!
        var f_2: Double!
        var f_3: Double!
        var f_4: Double!
        var f_5: Double!
        var f_6: Double!

        // estimated direction o the gyroscope error

        var J_11or24: Double!
        var J_12or23: Double!
        var J_13or22: Double!
        var J_14or21: Double!

        var J_32: Double!
        var J_33: Double!
        var J_41: Double!
        var J_42: Double!
        var J_43: Double!
        var J_44: Double!
        var J_51: Double!
        var J_52: Double!
        var J_53: Double!
        var J_54: Double!
        var J_61: Double!
        var J_62: Double!
        var J_63: Double!
        var J_64: Double!

        // estimated direction of the gyroscope error
        //  can't be initially unwrapped cause compiler tells warning:
        //      Binary operator '-=' cannot be applied to operands of type 'Double!' and 'Double'

        var SEqHatDot_1 = 0.0
        var SEqHatDot_2 = 0.0
        var SEqHatDot_3 = 0.0
        var SEqHatDot_4 = 0.0

        // estimated direction o the gyroscope error (angular)

        var w_err_x: Double!
        var w_err_y: Double!
        var w_err_z: Double!

        // computed flux in the earth frame

        var h_x: Double!
        var h_y: Double!
        var h_z: Double!

        // axulirary variables to avoid repeated calculations

        let halfSEq_1 = 0.5 * SEq_1
        let halfSEq_2 = 0.5 * SEq_2
        let halfSEq_3 = 0.5 * SEq_3
        let halfSEq_4 = 0.5 * SEq_4

        let twoSEq_1 = 2.0 * SEq_1
        let twoSEq_2 = 2.0 * SEq_2
        let twoSEq_3 = 2.0 * SEq_3
        let twoSEq_4 = 2.0 * SEq_4

        let twob_x = 2.0 * b_x
        let twob_z = 2.0 * b_z

        let twob_xSEq_1 = 2.0 * b_x * SEq_1
        let twob_xSEq_2 = 2.0 * b_x * SEq_2
        let twob_xSEq_3 = 2.0 * b_x * SEq_3
        let twob_xSEq_4 = 2.0 * b_x * SEq_4

        let twob_zSEq_1 = 2.0 * b_z * SEq_1
        let twob_zSEq_2 = 2.0 * b_z * SEq_2
        let twob_zSEq_3 = 2.0 * b_z * SEq_3
        let twob_zSEq_4 = 2.0 * b_z * SEq_4

        var SEq_1SEq_2: Double!
        var SEq_1SEq_3 = SEq_1 * SEq_3
        var SEq_1SEq_4: Double!
        var SEq_2SEq_3: Double!
        var SEq_2SEq_4 = SEq_2 * SEq_4
        var SEq_3SEq_4: Double!

        let twom_x = 2.0 * m_x
        let twom_y = 2.0 * m_y
        let twom_z = 2.0 * m_z

        // normalise the accelerometer measurement

        norm = sqrt(a_x * a_x + a_y * a_y + a_z * a_z)
        mutable_a_x /= norm
        mutable_a_y /= norm
        mutable_a_z /= norm

        // normalise the magnetometer measurement

        norm = sqrt(m_x * m_x + m_y * m_y + m_z * m_z)
        mutable_m_x /= norm
        mutable_m_y /= norm
        mutable_m_z /= norm

        // compute the objective function and Jacobian

        f_1 = twoSEq_2 * SEq_4 - twoSEq_1 * SEq_3 - a_x
        f_2 = twoSEq_1 * SEq_2 + twoSEq_3 * SEq_4 - a_y
        f_3 = 1.0 - twoSEq_2 * SEq_2 - twoSEq_3 * SEq_3 - a_z
        f_4 = twob_x * (0.5 - SEq_3 * SEq_3 - SEq_4 * SEq_4) + twob_z * (SEq_2SEq_4 - SEq_1SEq_3) - m_x
        f_5 = twob_x * (SEq_2 * SEq_3 - SEq_1 * SEq_4) + twob_z * (SEq_1 * SEq_2 + SEq_3 * SEq_4) - m_y
        f_6 = twob_x * (SEq_1SEq_3 + SEq_2SEq_4) + twob_z * (0.5 - SEq_2 * SEq_2 - SEq_3 * SEq_3) - m_z

        // J_11 negated in matrix multiplication

        J_11or24 = twoSEq_3
        J_12or23 = 2.0 * SEq_4

        // J_12 negated in matrix multiplication

        J_13or22 = twoSEq_1
        J_14or21 = twoSEq_2

        // negated in matrix multiplication

        J_32 = 2.0 * J_14or21
        J_33 = 2.0 * J_11or24

        J_41 = twob_zSEq_3
        J_42 = twob_zSEq_4
        J_43 = 2.0 * twob_xSEq_3 + twob_zSEq_1
        J_44 = 2.0 * twob_xSEq_4 - twob_zSEq_2

        J_51 = twob_xSEq_4 - twob_zSEq_2
        J_52 = twob_xSEq_3 + twob_zSEq_1
        J_53 = twob_xSEq_2 + twob_zSEq_4
        J_54 = twob_xSEq_1 - twob_zSEq_3

        J_61 = twob_xSEq_3
        J_62 = twob_xSEq_4 - 2.0 * twob_zSEq_2
        J_63 = twob_xSEq_1 - 2.0 * twob_zSEq_3
        J_64 = twob_xSEq_2

        // compute the gradient (matrix multiplication)

        // following expression solves following compiler warning:
        //  expression was too complex to be solved in reasonable time; consider breaking up the expression into distinct sub-expressions

        SEqHatDot_1 = J_14or21 * f_2
        SEqHatDot_1 -= J_11or24 * f_1
        SEqHatDot_1 -= J_41 * f_4
        SEqHatDot_1 -= J_51 * f_5
        SEqHatDot_1 += J_61 * f_6

        SEqHatDot_2 = J_12or23 * f_1
        SEqHatDot_2 += J_13or22 * f_2
        SEqHatDot_2 -= J_32 * f_3
        SEqHatDot_2 -= J_42 * f_4
        SEqHatDot_2 += J_52 * f_5
        SEqHatDot_2 += J_62 * f_6

        SEqHatDot_3 = J_12or23 * f_2
        SEqHatDot_3 -= J_33 * f_3
        SEqHatDot_3 -= J_13or22 * f_1
        SEqHatDot_3 -= J_43 * f_4
        SEqHatDot_3 += J_53 * f_5
        SEqHatDot_3 += J_63 * f_6

        SEqHatDot_4 = J_14or21 * f_1
        SEqHatDot_4 += J_11or24 * f_2
        SEqHatDot_4 -= J_44 * f_4
        SEqHatDot_4 -= J_54 * f_5
        SEqHatDot_4 += J_64 * f_6

        // normalise the gradient to estimate direction o the gyroscope error

        norm = sqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3 + SEqHatDot_4 * SEqHatDot_4)
        SEqHatDot_1 = SEqHatDot_1 / norm
        SEqHatDot_2 = SEqHatDot_2 / norm
        SEqHatDot_3 = SEqHatDot_3 / norm
        SEqHatDot_4 = SEqHatDot_4 / norm

        // compute angular estimated direction o the gyroscope error

        w_err_x = twoSEq_1 * SEqHatDot_2 - twoSEq_2 * SEqHatDot_1 - twoSEq_3 * SEqHatDot_4 + twoSEq_4 * SEqHatDot_3
        w_err_y = twoSEq_1 * SEqHatDot_3 + twoSEq_2 * SEqHatDot_4 - twoSEq_3 * SEqHatDot_1 - twoSEq_4 * SEqHatDot_2
        w_err_z = twoSEq_1 * SEqHatDot_4 - twoSEq_2 * SEqHatDot_3 + twoSEq_3 * SEqHatDot_2 - twoSEq_4 * SEqHatDot_1

        // compute and remove the gyroscope baises

        w_bx += w_err_x * deltat * zeta
        w_by += w_err_y * deltat * zeta
        w_bz += w_err_z * deltat * zeta
        mutable_w_x -= w_bx
        mutable_w_y -= w_by
        mutable_w_z -= w_bz

        // compute the quaternion rate measured by gyroscopes

        SEqDot_omega_1 = -halfSEq_2 * w_x - halfSEq_3 * w_y - halfSEq_4 * w_z
        SEqDot_omega_2 = halfSEq_1 * w_x + halfSEq_3 * w_z - halfSEq_4 * w_y
        SEqDot_omega_3 = halfSEq_1 * w_y - halfSEq_2 * w_z + halfSEq_4 * w_x
        SEqDot_omega_4 = halfSEq_1 * w_z + halfSEq_2 * w_y - halfSEq_3 * w_x

        // compute then integrate the estimated quaternion rate

        SEq_1 += (SEqDot_omega_1 - (beta * SEqHatDot_1)) * deltat
        SEq_2 += (SEqDot_omega_2 - (beta * SEqHatDot_2)) * deltat
        SEq_3 += (SEqDot_omega_3 - (beta * SEqHatDot_3)) * deltat
        SEq_4 += (SEqDot_omega_4 - (beta * SEqHatDot_4)) * deltat

        // normalise quaternion

        norm = sqrt(SEq_1 * SEq_1 + SEq_2 * SEq_2 + SEq_3 * SEq_3 + SEq_4 * SEq_4)
        SEq_1 /= norm
        SEq_2 /= norm
        SEq_3 /= norm
        SEq_4 /= norm

        // compute flux in the earth frame

        SEq_1SEq_2 = SEq_1 * SEq_2
        SEq_1SEq_3 = SEq_1 * SEq_3
        SEq_1SEq_4 = SEq_1 * SEq_4
        SEq_3SEq_4 = SEq_3 * SEq_4
        SEq_2SEq_3 = SEq_2 * SEq_3
        SEq_2SEq_4 = SEq_2 * SEq_4

        h_x = twom_x * (0.5 - SEq_3 * SEq_3 - SEq_4 * SEq_4) + twom_y * (SEq_2SEq_3 - SEq_1SEq_4) + twom_z * (SEq_2SEq_4 + SEq_1SEq_3)
        h_y = twom_x * (SEq_2SEq_3 + SEq_1SEq_4) + twom_y * (0.5 - SEq_2 * SEq_2 - SEq_4 * SEq_4) + twom_z * (SEq_3SEq_4 - SEq_1SEq_2)
        h_z = twom_x * (SEq_2SEq_4 - SEq_1SEq_3) + twom_y * (SEq_3SEq_4 + SEq_1SEq_2) + twom_z * (0.5 - SEq_2 * SEq_2 - SEq_3 * SEq_3)

        // normalise the flux vector to have only components in the x and z

        b_x = sqrt((h_x * h_x) + (h_y * h_y))
        b_z = h_z
    }
}

