# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import pandas as pd

class LogUtilities:
    def __init__(self, log_filename):
        with open(log_filename) as f:
            self.lines = f.readlines()

        log_data_type = ['time', 'dem_walltime']
        self.log_dataframe = pd.DataFrame(columns=log_data_type)
        self.log_dataframe['time'] = self.get_times()
        self.total_walltime = []

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

    def dem_performance_monitoring(self, n_char_0=49, n_char=8):
        walltime_keyword = "Total wallclock"
        walltime = [0]
        dof_time = True

        for line in self.lines:
            if walltime_keyword in line:
                if dof_time == True:
                    dof_time = False
                else:
                    if walltime[-1] > 999:

                        a = 1
                    walltime.append(float(line[n_char_0:n_char_0 + n_char].strip()))

        self.log_dataframe['dem_walltime'] = walltime
        self.total_walltime.append(walltime[-1])

        print("DEM time monitoring done.")

