import matplotlib.pyplot as plt

# File path
file_path = 'tst.out'

# Data containers
ntime = []
res_u, res_v, res_p, res_t = [], [], [], []

with open(file_path, 'r') as file:
    lines = file.readlines()[10:]  # Skip first 10 lines

    for line in lines:
        parts = line.split()
        if len(parts) >= 6:  # Must have at least 6 columns (NTIME ITER U V P ...)
            try:
                # Parse values first
                ntime_val = int(parts[0])         # NTIME
                resu_val = float(parts[2])        # Residual U
                resv_val = float(parts[3])        # Residual V
                resp_val = float(parts[4])        # Residual P
                rest_val = float(parts[5])        # Residual T

                # Only append if all parsing successful
                ntime.append(ntime_val)
                res_u.append(resu_val)
                res_v.append(resv_val)
                res_p.append(resp_val)
                res_t.append(rest_val)

            except ValueError:
                continue  # Skip malformed lines

# Confirm everything is same length
print(f"Lengths -> NTIME: {len(ntime)}, ResU: {len(res_u)}, ResV: {len(res_v)}, ResP: {len(res_p)}, ResT: {len(res_t)}")
print(ntime)
# Plotting
plt.figure(figsize=(10, 6))
plt.plot(ntime, res_u, label='Residual U')
plt.plot(ntime, res_v, label='Residual V')
plt.plot(ntime, res_p, label='Residual P')
plt.plot(ntime, res_t, label='Residual T')
plt.yscale('log')
plt.xlabel('Time Step (NTIME)')
plt.ylabel('Residual')
plt.title('Residuals vs Time Step')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
