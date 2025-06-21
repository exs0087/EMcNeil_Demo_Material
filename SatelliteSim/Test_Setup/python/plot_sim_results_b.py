#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main(csv_path, out_dir):
    # load
    df = pd.read_csv(csv_path)
    t  = df['time'].to_numpy()

    # extract columns
    w1 = df['wx'].to_numpy()
    w2 = df['wy'].to_numpy()
    w3 = df['wz'].to_numpy()
    q1 = df['q1'].to_numpy()
    q2 = df['q2'].to_numpy()
    q3 = df['q3'].to_numpy()
    q4 = df['q4'].to_numpy()
    m1 = df['mx'].to_numpy()
    m2 = df['my'].to_numpy()
    m3 = df['mz'].to_numpy()

    # magnitude of angular velocity
    w_mag = np.sqrt(w1**2 + w2**2 + w3**2)

    # 1) Quaternion
    plt.figure()
    plt.plot(t, q1, '-', label='q1')
    plt.plot(t, q2, '-', label='q2')
    plt.plot(t, q3, '-', label='q3')
    plt.plot(t, q4, '-', label='q4')
    plt.xlim(0, t[-1])
    plt.ylim(-1, 1)
    plt.grid(True)
    plt.legend()
    plt.title('Erin McNeil: HW 3, Prob 3,b: Quaternion Plot')
    plt.xlabel('time (s)')
    plt.ylabel('quaternion values')
    plt.savefig(f"{out_dir}/HW3_p3b_quaternion.png")

    # 2) Angular velocity
    plt.figure()
    plt.plot(t, w1, 'r-*', markevery=500, label='w1')
    plt.plot(t, w2, 'b-',        label='w2')
    plt.plot(t, w3, 'm-',        label='w3')
    plt.xlim(0, t[-1])
    plt.ylim(-0.4,  0.4)
    plt.grid(True)
    plt.legend()
    plt.title('Erin McNeil: HW 3, Prob 3,b: Angular Velocity Plot')
    plt.xlabel('time (s)')
    plt.ylabel('angular velocity (rad/sec)')
    plt.savefig(f"{out_dir}/HW3_p3b_angvel.png")

    # 3) Magnitude of angular velocity
    plt.figure()
    plt.plot(t, w_mag, 'm-')
    plt.xlim(0, t[-1])
    plt.grid(True)
    plt.legend(['magnitude of angular velocity'])
    plt.title('Erin McNeil: HW 3, Prob 3,b: Magnitude of Angular Velocity Plot')
    plt.xlabel('time (s)')
    plt.ylabel('magnitude of angular velocity (rad/sec)')
    plt.savefig(f"{out_dir}/HW3_p3b_angvel_mag.png")

    # 4) Magnetic dipole moment
    plt.figure()
    plt.plot(t, m1, 'r-*', markevery=500, label='mb1')
    plt.plot(t, m2, 'b-',        label='mb2')
    plt.plot(t, m3, 'm-',        label='mb3')
    plt.xlim(0, t[-1])
    plt.grid(True)
    plt.legend()
    plt.title('Erin McNeil: HW 3, Prob 3,b: Magnetic Dipole Moment Plot')
    plt.xlabel('time (s)')
    plt.ylabel('Magnetic Dipole Moment (A·m²)')
    plt.savefig(f"{out_dir}/HW3_p3b_dipoles.png")

    print("All HW3-b plots saved to", out_dir)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Plot HW3 Prob 3b results")
    p.add_argument("csv",
                   nargs="?",
                   default=r"C:\Users\ees00\Demo_Github\SatelliteSim\Test_Setup\cpp\sim_results.csv",
                   help="path to sim_results.csv")
    p.add_argument("-o","--out",
                   default="../cpp/plots",
                   help="output directory for PNGs")
    args = p.parse_args()
    main(args.csv, args.out)
