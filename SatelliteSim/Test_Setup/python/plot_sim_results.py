#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main(csv_path, out_dir):
    # 1) Load the data
    df = pd.read_csv(csv_path)
    t = df['time']

    # 2) Plot angular rates
    plt.figure()
    plt.plot(t, df['wx'], label='ωₓ')
    plt.plot(t, df['wy'], label='ω_y')
    plt.plot(t, df['wz'], label='ω_z')
    plt.xlabel('Time [s]')
    plt.ylabel('Angular rate [rad/s]')
    plt.legend()
    plt.title('Body‐rates vs. Time')
    plt.grid(True)
    plt.savefig(f"{out_dir}/body_rates.png")

    # 3) Plot quaternion components
    plt.figure()
    for q in ['q1','q2','q3','q4']:
        plt.plot(t, df[q], label=q)
    plt.xlabel('Time [s]')
    plt.ylabel('Quaternion component')
    plt.legend()
    plt.title('Attitude Quaternion vs. Time')
    plt.grid(True)
    plt.savefig(f"{out_dir}/quaternion.png")

    # 4) Plot magnetic field in body frame
    plt.figure()
    plt.plot(t, df['Bx'], label='Bₓ')
    plt.plot(t, df['By'], label='B_y')
    plt.plot(t, df['Bz'], label='B_z')
    plt.xlabel('Time [s]')
    plt.ylabel('Magnetic field [T]')
    plt.legend()
    plt.title('Magnetic Field vs. Time')
    plt.grid(True)
    plt.savefig(f"{out_dir}/mag_field.png")

    # 5) Plot commanded dipoles
    plt.figure()
    plt.plot(t, df['mx'], label='mₓ')
    plt.plot(t, df['my'], label='m_y')
    plt.plot(t, df['mz'], label='m_z')
    plt.xlabel('Time [s]')
    plt.ylabel('Dipole moment [A·m²]')
    plt.legend()
    plt.title('Commanded Magnetic Dipoles vs. Time')
    plt.grid(True)
    plt.savefig(f"{out_dir}/dipoles.png")

    print("Plots saved to", out_dir)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Plot C++ sim results")
    p.add_argument("csv",
                   nargs="?",
                   # raw string—no escape sequences processed
                   default=r"C:\Users\ees00\Demo_Github\SatelliteSim\Test_Setup\cpp\sim_results.csv",
    help="path to sim_results.csv")
    p.add_argument("-o", "--out", default=".", help="Output directory for PNGs")
    args = p.parse_args()
    main(args.csv, args.out)
