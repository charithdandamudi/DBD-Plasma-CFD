import matplotlib.pyplot as plt
import matplotlib.animation as animation

file_path = 'tst.out'

# Create figure and axes
fig, ax = plt.subplots(figsize=(10, 7))

# Create white space above main plot
plt.subplots_adjust(top=0.87)  # adjust this to control the white space height

# Residual plot setup
line_u, = ax.plot([], [], label='Residual U')
line_v, = ax.plot([], [], label='Residual V')
line_p, = ax.plot([], [], label='Residual P')
line_t, = ax.plot([], [], label='Residual T')

ax.set_yscale('log')
ax.set_xlabel('Time Step (NTIME)')
ax.set_ylabel('Residual')
ax.set_title('Live Residuals vs Time Step')
ax.legend()
ax.grid(True)

# Get the right and top edges of the plot
plot_box = ax.get_position()
plot_right = plot_box.x1
plot_top = plot_box.y1

# Layout parameters
bar_height = 0.015
bar_y = plot_top + 0.07
text_width = 0.07    # adjust this if text is wider
spacing = 0.01       # gap between bar and text
bar_width = 0.13     # width of the progress bar
text_x = plot_right
bar_x = text_x - text_width - spacing - bar_width

# Set up bar axes
pax = fig.add_axes([bar_x, bar_y, bar_width, bar_height])
pax.set_xlim(0, 1)
pax.set_ylim(0, 1)
pax.axis('off')
pax.barh(0.5, width=1.0, height=0.4, color='#E0E0E0', zorder=0)
progress_bar = pax.barh(0.5, width=0.0, height=0.4, color='#4A90E2', zorder=1)[0]

# Set up anchored text
progress_text = fig.text(text_x, bar_y + bar_height / 2,
                         '', ha='right', va='center',
                         fontsize=9, color='#333333',
                         transform=fig.transFigure)



# Data containers
ntime_limit = None
freeze_updates = False
last_seen_nt = 0
ntime_vals, res_u, res_v, res_p, res_t = [], [], [], [], []

def read_ntime_limit():
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            return int(lines[1].split()[0])
    except Exception:
        return None

def init():
    ax.set_xlim(0, 10)
    ax.set_ylim(1e-10, 1e0)
    progress_bar.set_width(0)
    progress_text.set_text('')
    return line_u, line_v, line_p, line_t, progress_bar, progress_text

def update(frame):
    global ntime_limit, freeze_updates, last_seen_nt
    global ntime_vals, res_u, res_v, res_p, res_t

    if freeze_updates:
        return line_u, line_v, line_p, line_t, progress_bar, progress_text

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        if ntime_limit is None and len(lines) > 1:
            ntime_limit = int(lines[1].split()[0])

        lines = lines[10:]  # Skip header

        # Early check: skip update if NT < 2 only
        nt_values = [int(line.split()[0]) for line in lines if line.strip() and line.split()[0].isdigit()]
        if not any(nt >= 2 for nt in nt_values):
            # print("Waiting for NT >= 2...")
            return line_u, line_v, line_p, line_t, progress_bar, progress_text

        new_ntime, new_u, new_v, new_p, new_t = [], [], [], [], []

        for line in lines:
            parts = line.split()
            if len(parts) >= 6:
                try:
                    nt = int(parts[0])
                    new_ntime.append(nt)
                    new_u.append(float(parts[2]))
                    new_v.append(float(parts[3]))
                    new_p.append(float(parts[4]))
                    new_t.append(float(parts[5]))

                    if nt >= ntime_limit:
                        freeze_updates = True
                        print(f"\nReached NTIME = {ntime_limit}. Freezing Plot Updates.")
                        break

                except ValueError:
                    continue

        if not new_ntime:
            return line_u, line_v, line_p, line_t, progress_bar, progress_text

        current_nt = new_ntime[-1]
        if current_nt != last_seen_nt:
            last_seen_nt = current_nt

            ntime_vals = new_ntime
            res_u = new_u
            res_v = new_v
            res_p = new_p
            res_t = new_t

            line_u.set_data(ntime_vals, res_u)
            line_v.set_data(ntime_vals, res_v)
            line_p.set_data(ntime_vals, res_p)
            line_t.set_data(ntime_vals, res_t)

            ax.set_xlim(0, max(ntime_vals) + 1)
            ax.set_ylim(1e-10, max(max(res_u + res_v + res_p + res_t), 1e-8))

            # Compute progress percentage
            progress = current_nt / ntime_limit
            progress_bar.set_width(progress)

            # Change color when progress hits 100%
            if progress >= 1.0:
                progress_bar.set_color('green')
            else:
                progress_bar.set_color('#4A90E2')  # original blue


            # Format progress string
            progress_str = f"{int(progress * 100)}%  ({current_nt}/{ntime_limit})"
            progress_text.set_text(progress_str)

            # Estimate text width
            text_len = len(progress_str)
            text_width = text_len * 0.0075

            # Fixed total width for bar + spacing + text
            total_width = 0.2
            spacing = 0.008

            # Shrink bar width based on text size
            bar_width = max(total_width - text_width - spacing, 0.05)  # clamp to avoid disappearing bar
            bar_x = ax.get_position().x1 - text_width - spacing - bar_width

            # Reposition bar and set width
            pax.set_position([bar_x, pax.get_position().y0, bar_width, pax.get_position().height])


    except Exception:
        pass

    return line_u, line_v, line_p, line_t, progress_bar, progress_text

ani = animation.FuncAnimation(fig, update, init_func=init, interval=1000, blit=False)
plt.show()
