# Sort all data given reference array
def sort_by_array(self, reference_array_name="ID"):
    """
    Sorts all self.df according to a reference array:

    Parameter:

    :param reference_array_name = "ID" -> String with name of reference array.
    "ID" is used as default for particles, but any other 1D array can be
    used.
    """

    if self.sorted:
        return

    global sort_by_array_loop

    if self.df_available:
        def sort_by_array_loop(i):
            self.df[i].points = self.df[i].points[self.df[i][reference_array_name].argsort()]
            for name in self.df[0].array_names:
                self.df[i][name] = self.df[i][name][self.df[i][reference_array_name].argsort()]

    else:
        def sort_by_array_loop(i):
            df = self.get_df(i)
            df.points = df.points[df[reference_array_name].argsort()]
            for name in df.array_names:
                df[name] = df[name][df[reference_array_name].argsort()]

            df.save(f'{self.path_output}/{self.list_vtu[i]}')

    self.parallel_run(sort_by_array_loop, range(len(self.list_vtu)),
                      tqdm_desc=f"Sorting dataframe by {reference_array_name}")

    self.sorted = True
