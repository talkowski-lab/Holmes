import numpy as np


class Event(object):
    def __init__(self, data):
        self.name, self.sample, self.var_type = data[0:3]
        self.qualA, self.qualB, self.qual_res = [float(x) for x in data[3:6]]
        self.uniqA, self.uniqB, self.uniq_res = [float(x) for x in data[6:9]]
        self.spanA, self.spanB, self.span_res = [float(x) for x in data[9:12]]
        self.global_cov = float(data[12])
        self.localA, self.localB, self.local_res = [float(x) for x in data[13:16]]
        self.gcA, self.gcB, self.gc_res = [float(x) for x in data[16:19]]
        self.alignA, self.alignB, self.align_res = [float(x) for x in data[19:22]]

    @classmethod
    def init_link(cls, clusters, name):
        qualA = np.mean([x.qualA for x in clusters])
        qualB = np.mean([x.qualB for x in clusters])
        qual_res = np.mean([x.qual_res for x in clusters])

        uniqA = np.mean([x.uniqA for x in clusters])
        uniqB = np.mean([x.uniqB for x in clusters])
        uniq_res = np.mean([x.uniq_res for x in clusters])

        spanA = np.mean([abs(x.spanA) for x in clusters])
        spanB = np.mean([abs(x.spanB) for x in clusters])
        span_res = np.mean([x.span_res for x in clusters])

        global_cov = np.mean([x.global_cov for x in clusters])

        localA = np.mean([x.localA for x in clusters])
        localB = np.mean([x.localB for x in clusters])
        local_res = np.mean([x.local_res for x in clusters])

        gcA = np.mean([x.gcA for x in clusters])
        gcB = np.mean([x.gcB for x in clusters])
        gc_res = np.mean([x.gc_res for x in clusters])

        alignA = np.mean([x.alignA for x in clusters])
        alignB = np.mean([x.alignB for x in clusters])
        align_res = np.mean([x.align_res for x in clusters])

        var_types = [x.name.split('_')[-2] for x in clusters]

        # Determine variant type
        # If more than one type of variant cluster in event, set as 'complex'
        if len(set(var_types)) == 1:
            var_type = var_types[0]
            sample = clusters[0].name.split(var_type)[0].rstrip('_')
        else:
            var_type = 'complex'
            sample = clusters[0].name.split(var_types[0])[0].rstrip('_')

        return Event([name, sample, var_type,
                     qualA, qualB, qual_res,
                     uniqA, uniqB, uniq_res,
                     spanA, spanB, span_res,
                     global_cov,
                     localA, localB, local_res,
                     gcA, gcB, gc_res,
                     alignA, alignB, align_res])

    def __repr__(self):
        return ("%s\t%s\t%s\t"
                "%.2f\t%.2f\t%.2f\t%.2f\t"
                "%.2f\t%.2f\t"
                "%.2f\t%.2f\t%.2f\t"
                "%.2f\t"
                "%.2f\t%.2f\t%.2f\t"
                "%.2f\t%.2f\t%.2f\t"
                "%.2f\t%.2f\t%.2f\t" %
                (self.name, self.sample, self.var_type,
                self.qualA, self.qualB, self.uniqA, self.uniqB,
                self.qual_res, self.uniq_res,
                self.spanA, self.spanB, self.span_res,
                self.global_cov,
                self.localA, self.localB, self.local_res,
                self.gcA, self.gcB, self.gc_res,
                self.alignA, self.alignB, self.align_res))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
