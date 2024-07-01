import numpy as np
import pandas as pd

class LogUtilities:
    def __init__(self, log_filename):
        with open(log_filename) as f:
            self.lines = f.readlines()

        log_data_type = ['time', 'beta', 'fluid_velocity']
        self.log_dataframe = pd.DataFrame(columns=log_data_type)
        self.log_dataframe['time'] = self.get_times()

    def get_times(self):
        # Find the time in the log file through the keyword "Transient" and some hardcoded indices
        time_keyword = "Transient"
        time = np.array([0])

        for line in self.lines:
            if time_keyword in line:
                time_str = line[36:36 + 6].strip()
                time = np.append(time, float(time_str))
        return time

    def get_log_data(self):
        return self.log_dataframe

    def flow_monitoring(self):
        velocity_keyword = "Fluid space-average"
        beta_keyword = "Fluid beta force"

        average_velocity = np.array([0])
        beta = np.array([0])

        for line in self.lines:
            if velocity_keyword in line:
                average_velocity = np.append(average_velocity, float(line.split()[-1]))
            if beta_keyword in line:
                beta = np.append(beta, float(line.split()[-1]))

        self.log_dataframe['fluid_velocity'] = average_velocity
        self.log_dataframe['beta'] = beta

        print("Flow monitoring done.")