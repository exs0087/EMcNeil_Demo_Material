import numpy as np

def test_columns(df):
    assert list(df.columns) == [
        "time","wx","wy","wz",
        "q1","q2","q3","q4",
        "Bx","By","Bz",
        "mx","my","mz"
    ]

def test_length(df):
    assert len(df) == 16000 + 1

def test_initial_state(df):
    row0 = df.iloc[0]
    assert row0["time"] == 0
    assert (row0[["wx","wy","wz"]] == 0).all()
    assert (row0[["q1","q2","q3"]] == 0).all()
    assert row0["q4"] == 1

def test_quaternion_normalization(df):
    norms = np.sqrt(
        df["q1"]**2 + df["q2"]**2 + df["q3"]**2 + df["q4"]**2
    )
    assert np.allclose(norms, 1.0, atol=1e-6)

def test_B_field_magnitude(df):
    Bmag = np.sqrt(df["Bx"]**2 + df["By"]**2 + df["Bz"]**2)
    # all finite
    assert np.isfinite(Bmag).all()
    # within expected range
    assert (Bmag > 1e-9).all() and (Bmag < 1e-4).all()

def test_zero_dipole(df):
    assert (df[["mx","my","mz"]].values == 0).all()
