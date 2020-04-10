from pyarrow import parquet as pq
import numpy
import pandas

def _read(file, columns=None, skip_individuals=False, to_pandas=False, specific_individuals=None):
    if columns is None:
        columns = file.schema.names
    if not skip_individuals:
        columns = ["individual"]+[x for x in columns]
    if skip_individuals and specific_individuals is not None:
            raise RuntimeError("Unsupported combination")
    v = file.read(columns=columns)
    if to_pandas:
        v = v.to_pandas()
        if specific_individuals is not None:
            indexes = { str(i) for i in specific_individuals}
            v = v.loc[v.individual.isin(indexes)]
    else:
        if specific_individuals:
            mask = _individual_mask(v.column(0).to_pylist(), specific_individuals)

        v = {c.name:(numpy.array(c.to_pylist(), dtype=numpy.float32) if c.name != "individual" else numpy.array(c.to_pylist(), dtype=numpy.str)) for c in v}
        if specific_individuals:
            v = {k:d[mask] for k,d in v.items()}
    return v

def _individual_mask(individuals, specific_individuals):
    if specific_individuals:
        return [individuals.index(x) for x in specific_individuals]
    else:
        return None

class MultiFileGenoHandler:
    """
    This class is for loading parquet metadata and genotype files. Most of its
    functionality is meant to assist in the case that both metadata and genotype
    are split into 22 files, but it should be robust to the case where there is
    only one genotype or metadata file.
    """
    def __init__(self, features, metadata):
        """
        If either argument is a pattern for multiple files, it must be
        formattable with the argument 'chr'

        :param features: filepath (or pattern for multiple) for genotype files
        :param metadata: filepath (or pattern for multiple) for geno metadata
        """
        self.m_features = self.check_if_formattable(features)
        if self.m_features:
            self.features = self.format_chrom_file_names(features)
        else:
            self.features = [features]

        self.m_metadata = self.check_if_formattable(metadata)
        if self.m_metadata:
            self.metadata = self.format_chrom_file_names(metadata)
        else:
            self.metadata = [metadata]

    @staticmethod
    def format_chrom_file_names(s):
        l = [s.format(chr=i) for i in range(1, 23)]
        return l
    @staticmethod
    def check_if_formattable(s):
        matches = re.findall('{(.*?)}', s)
        if len(matches) > 0 and matches[0] == 'chr':
            return True
        else:
            return False

    def load_metadata(self, whitelist=None):
        df_lst = []
        for i in self.metadata:
            df_i = pq.read_table(i).to_pandas()
            if whitelist is not None:
                df_i = df_i[df_i.id.isin(whitelist)]
            df_lst.append(df_i)
        return pandas.concat(df_lst)

    def load_features(self, metadata, individuals, pandas=False):
        """
        :param metadata: pandas DataFrame with columns 'variant_id' and
                'chromosome'
        :param individuals: list. Individual IDs
        :param pandas: bool. default False. Whether the returned obj should be
                a pandas DataFrame
        :return: dict.
        """
        if self.m_features:
            return self._load_features_multiple(metadata, individuals, pandas)
        else:
            return self._load_features_single(metadata, individuals, pandas)

    def _load_features_single(self, metadata, individuals, pandas):
        dd =  _read(pq.ParquetFile(self.features[0]),
                    columns=[x for x in metadata.id],
                    specific_individuals=individuals,
                    to_pandas=pandas)
        logging.log(5, "Loaded {} features".format(len(dd) - 1))
        return dd

    def _load_features_multiple(self, metadata, individuals, pandas):
        df_lst = []
        i_ = list(individuals)
        for chr, group in metadata.groupby('chromosome'):
            chr_fp = self.features[chr - 1]
            chr_vars = list(group.id)
            chr_features = _read(pq.ParquetFile(chr_fp), chr_vars,
                                         specific_individuals=i_,
                                         to_pandas = True)
            df_lst.append(chr_features.set_index('individual'))
        while len(df_lst) > 1:
            df_lst[0].join(df_lst.pop(), how='inner')
        logging.log(5, "Loaded {} features".format(df_lst[0].shape[1]))
        if pandas:
            return df_lst[0]
        else:
            return df_lst[0].reset_index().to_dict(orient='list')
