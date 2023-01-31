import pandas as pd
import numpy as np
import Charts as Ch


class Plane:
    """Треугольный КЭ пластины полоской деформации."""

    def __init__(self, type, t, v, E):
        """Constructor."""
        self.type = type
        self.t = t
        self.v = v
        self.E = E
        self.Pl = pd.read_csv('PlaneData.csv', sep=';')
        self.npl = self.Pl.shape[0]

    def setdata(self, p):
        """Данные КЭ."""

        self.pldat = self.Pl.loc[p, :]
        self.gl = self.pldat[:3].to_numpy().astype(int)  # Узлы

    def prop(self, Ji, Jj, Jk):
        """Свойства элемента."""

        self.jn = np.array([[Ji.jntdat['X'], Ji.jntdat['Y']],
                            [Jj.jntdat['X'], Jj.jntdat['Y']],
                            [Jk.jntdat['X'], Jk.jntdat['Y']]])  # Координаты узлов элемента

        self.A = 1 / 2 * np.linalg.det(np.array([
            [1, self.jn[0, 0], self.jn[0, 1]],
            [1, self.jn[1, 0], self.jn[1, 1]],
            [1, self.jn[2, 0], self.jn[2, 1]]]))  # Площадь КЭ

    def mD(self, v, E):
        """Матрица коэффициентов закона Гука."""

        if self.type == 'strain':
            self.txt = 'Плоская деформация'
            self.D = E / ((1 - 2 * v) * (1 + v)) * \
                     np.matrix([[1 - v, v, 0],
                                [v, 1 - v, 0],
                                [0, 0, (1 - 2 * v) / 2]])
        else:
            self.txt = 'Плоское напряжение'
            self.D = E / (1 - v ** 2) * \
                     np.matrix([[1, v, 0],
                                [v, 1, 0],
                                [0, 0, (1 - v) / 2]])

    def bc(self, l):
        if l == 0:
            b = self.jn[1, 1] - self.jn[2, 1]
            c = self.jn[2, 0] - self.jn[1, 0]
        elif l == 1:
            b = self.jn[2, 1] - self.jn[0, 1]
            c = self.jn[0, 0] - self.jn[2, 0]
        else:
            b = self.jn[0, 1] - self.jn[1, 1]
            c = self.jn[1, 0] - self.jn[0, 0]
        return [b, c]

    def mBr(self):
        """Матрица градиентов КЭ пластины."""

        self.Br = 1 / (2 * self.A) * np.matrix(
            [[self.bc(0)[0], 0, self.bc(1)[0], 0, self.bc(2)[0], 0],
             [0, self.bc(0)[1], 0, self.bc(1)[1], 0, self.bc(2)[1]],
             [self.bc(0)[1], self.bc(0)[0], self.bc(1)[1], self.bc(1)[0], self.bc(2)[1], self.bc(2)[0]]])

    def mk(self):
        """Матрица жесткости КЭ пластины."""

        Jointi.setdata(self.pldat['i'])
        Jointj.setdata(self.pldat['j'])
        Jointk.setdata(self.pldat['k'])
        Planei.prop(Jointi, Jointj, Jointk)
        Planei.mD(self.v, self.E)
        Planei.mBr()
        self.kr = (np.transpose(Planei.Br) * Planei.D * Planei.Br) * self.A * self.t


class Joint:
    """Узел КЭ пластины."""

    def __init__(self):
        """Constructor."""

        self.Jn = pd.read_csv('Joints.csv', sep=';')
        self.vert = self.Jn.values[:, :2]
        self.nj = self.Jn.shape[0]

    def setdata(self, j):
        """Данные узла КЭ."""

        self.jntdat = self.Jn.loc[j, :]


class Model:
    """Расчетная модель."""

    def mK(self):
        """Матрица жесткости в глобальной СК."""

        self.jc = Jointi.Jn.values[:, :2]  # Координаты узлов
        self.K = np.zeros((Jointi.nj * 2, Jointi.nj * 2))
        self.xy = np.zeros((3, 2))
        for i in range(0, Planei.npl):
            Planei.setdata(i)
            Planei.mk()
            self.K = self.K + self.assembly()
            self.xy = np.append(self.xy, Planei.jn, axis=0)
        self.xy = self.xy[3:]
        self.jntloads()

    def assembly(self):
        """Сборка глобальной матрицы жесткости."""

        Ki = np.zeros((Jointi.nj * 2, Jointi.nj * 2))
        ii = 2 * Planei.gl
        jj = ii + 1
        ij = [ii[0], jj[0], ii[1], jj[1], ii[2], jj[2]]
        for i in range(6):
            for j in range(6):
                Ki[ij[i], ij[j]] += Planei.kr[i, j]
        return Ki

    def jntloads(self):
        """Узловые нагрузки."""

        self.Fp = np.zeros((Jointi.nj * 2, 1))
        for i in range(0, Jointi.nj):
            self.Fp[i * 2] = Jointi.Jn['Fx'][i]
            self.Fp[i * 2 + 1] = Jointi.Jn['Fy'][i]

    def rstrs(self):
        """Граничные условия."""

        for i in range(0, Jointi.nj):
            rstr = Jointi.Jn[['UX', 'UY']].iloc[i]
            if sum(rstr) != 0:
                for j in range(0, 2):
                    if rstr[j] != 0:
                        rs = i * 2 + j
                        for k in range(0, Jointi.nj * 2):
                            self.K[rs, k] = 0
                            self.K[k, rs] = 0
                        self.K[rs, rs] = 1
                        self.Fp[rs] = 0
            if rstr[0] == 1 or rstr[1] == 1 and rstr[2] == 0:
                self.K[i * 3 + 2, i * 3 + 2] = 1

    def sol(self):
        """Решение."""

        self.U = np.dot(np.linalg.inv(self.K), self.Fp)


if __name__ == "__main__":
    type = 'strain'
    t = 0.1
    v = 0.3
    E = 100000
    Planei = Plane(type, t, v, E)
    Jointi = Joint()
    Jointj = Joint()
    Jointk = Joint()
    Model1 = Model()
    Model1.mK()
    Model1.rstrs()
    Model1.sol()
    pd.DataFrame(Model1.K).to_csv('K.csv', sep=';')
    print('Перемещения узлов')
    print(pd.DataFrame(np.round(Model1.U.reshape(Jointi.nj, 2) * 1000, 6), columns=['UX[мм]', 'UY[мм]']))
    pd.DataFrame(Model1.U.reshape(Jointi.nj, 2), columns=['UX', 'UY']).to_csv('U.csv', sep=';')
    Ch.geom(Model1.xy, Jointi.vert, Planei.npl, Jointi.nj, 'Элементы')
    k = 10  # Коэффициент масштабирования перемещений
    dc = Model1.jc + Model1.U.reshape(Jointi.nj, 2) * k
    xyd = np.zeros((3, 2))
    for i in range(0, Planei.npl):
        Planei.setdata(i)
        gl = Planei.gl
        joints = np.array([[dc[gl[0], 0], dc[gl[0], 1]],  # Координаты узлов элемента
                           [dc[gl[1], 0], dc[gl[1], 1]],
                           [dc[gl[2], 0], dc[gl[2], 1]]])
        xyd = np.append(xyd, joints, axis=0)
    xyd = xyd[3:]  # Координаты узлов элемента + перемещения
    Ch.geom(xyd, dc, Planei.npl, Jointi.nj, Planei.txt)
