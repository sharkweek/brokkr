from pandas import DataFrame
from lasso.dyna import Binout

class Extended_Binout(Binout):
    def read(self, *args):
        """Extends the `read` functionality of `Binout`"""
        if args[0:2] == ('matsum', 'legend'):
            return self._matsum_legend_df()
        else:
            super.read(*args)

    def _matsum_legend_df(self):
        """reformat the `matsum` legend as a pandas `DataFrame`"""
        legend = super().read('matsum', 'legend')
        titles = [legend[i:i + 80].strip() for i in range(0, len(legend), 80)]
        parts = [list(x) for x
                 in zip(range(0, len(titles)),
                        super().read('matsum', 'ids'),
                        titles)]
        df = pd.DataFrame({
            'pid': super().read('matsum', 'ids'),
            'title': [legend[i:i + 80].strip()
                      for i in range(0, len(legend), 80)]
        })

        return df
