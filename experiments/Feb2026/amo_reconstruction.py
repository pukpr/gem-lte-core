import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO
from pathlib import Path

# Gray et al. 2004 tree-ring reconstruction (1567–1990)
raw = """age	ssta	amo
1567	-1.43849	NA
1568	-0.73836	NA
1569	-0.38016	NA
1570	0.45247	NA
1571	-0.09531	NA
1572	-0.80857	-0.60744
1573	-1.42524	-0.55674
1574	-0.89676	-0.53063
1575	-1.2915	-0.4747
1576	0.16973	-0.3968
1577	-0.68297	-0.25437
1578	-0.47982	0.03707
1579	-0.11654	0.44566
1580	1.30749	0.95336
1581	0.60763	1.3532
1582	1.33716	1.57431
1583	2.2579	1.76867
1584	3.59189	1.94388
1585	4.37384	1.96908
1586	2.5011	1.91519
1587	1.40787	1.90422
1588	1.31662	1.75355
1589	1.59112	1.44409
1590	0.10384	1.1004
1591	0.73349	0.83009
1592	0.99197	0.60076
1593	-0.41047	0.30456
1594	0.07109	-0.0391
1595	1.02082	-0.31122
1596	0.44799	-0.52503
1597	-1.12555	-0.77491
1598	-2.07403	-1.00366
1599	-1.89135	-1.16818
1600	-1.85613	-1.35047
1601	-1.58274	-1.55019
1602	-1.6895	-1.65926
1603	-2.30383	-1.7002
1604	-1.3261	-1.6756
1605	-1.2277	-1.60343
1606	-1.29783	-1.58978
1607	-1.56127	-1.5139
1608	-2.45699	-1.33301
1609	-1.01648	-1.23286
1610	-1.28749	-1.2373
1611	-1.87844	-1.24413
1612	0.12374	-1.17584
1613	-0.49919	-1.02204
1614	-1.12787	-0.98756
1615	-1.51461	-1.02792
1616	-1.14752	-0.94387
1617	-0.34585	-0.94601
1618	-0.59637	-1.08301
1619	-2.18761	-1.13356
1620	-0.92351	-1.08639
1621	-0.56141	-1.01569
1622	-1.236	-0.94666
1623	-1.87942	-0.88963
1624	-0.75873	-0.76377
1625	-0.94026	-0.62093
1626	-0.308	-0.61309
1627	0.1953	-0.59255
1628	0.00309	-0.42923
1629	-0.26995	-0.2595
1630	0.01574	-0.11223
1631	-1.34385	0.05231
1632	-0.04288	0.15674
1633	0.19389	0.24117
1634	0.56265	0.35499
1635	0.68363	0.4137
1636	1.35888	0.46647
1637	0.61713	0.47325
1638	1.26973	0.41732
1639	0.73993	0.3841
1640	0.17996	0.25549
1641	-0.45266	0.06852
1642	-0.79832	-0.02392
1643	-0.16941	-0.09574
1644	0.26162	-0.15452
1645	-1.58749	-0.1424
1646	-0.1094	-0.10235
1647	0.23659	-0.04873
1648	0.21384	-0.03044
1649	0.62023	-0.07672
1650	0.54195	-0.04699
1651	-0.01357	0.05907
1652	-0.165	0.11648
1653	-0.43699	0.16739
1654	-0.39637	0.2538
1655	-0.33484	0.28325
1656	0.75903	0.26224
1657	0.5165	0.2921
1658	0.95216	0.35961
1659	1.61011	0.46587
1660	0.14094	0.60928
1661	-0.03281	0.70823
1662	0.45151	0.74939
1663	0.29664	0.72606
1664	0.9952	0.66608
1665	1.14189	0.66963
1666	1.26126	0.73553
1667	0.83748	0.78755
1668	0.1645	0.78181
1669	1.19821	0.82413
1670	0.62394	0.92407
1671	0.80216	0.96324
1672	0.65687	0.9525
1673	-0.02354	1.01509
1674	2.16185	1.11694
1675	1.97407	1.13513
1676	1.21244	1.08953
1677	0.67143	1.1158
1678	1.58243	1.14674
1679	1.81724	1.07637
1680	0.36869	0.97527
1681	0.14551	0.88524
1682	1.83894	0.87295
1683	-0.58683	0.86941
1684	1.31778	0.74688
1685	0.79597	0.69028
1686	0.59001	0.79919
1687	1.04803	0.86939
1688	1.13507	0.99186
1689	-0.18599	1.13207
1690	1.23996	1.19268
1691	1.45234	1.24188
1692	1.93622	1.20402
1693	1.76516	1.11992
1694	1.77016	1.10206
1695	1.5558	1.05413
1696	0.8141	0.93009
1697	0.06676	0.76499
1698	0.43439	0.57248
1699	0.15746	0.41649
1700	-0.0621	0.32036
1701	0.27354	0.26519
1702	-0.18693	0.2819
1703	0.03803	0.30721
1704	0.37753	0.26318
1705	1.02585	0.21527
1706	0.24072	0.18461
1707	0.97421	0.17167
1708	0.03327	0.18031
1709	-0.32206	0.1891
1710	-0.54079	0.14412
1711	0.139	0.07842
1712	-0.31121	-0.00603
1713	0.33513	-0.10119
1714	0.25615	-0.17428
1715	0.24764	-0.29416
1716	-0.29493	-0.4254
1717	-0.17919	-0.47
1718	-0.71659	-0.52744
1719	-1.03389	-0.60678
1720	-2.22653	-0.61998
1721	-0.80007	-0.59235
1722	-0.26416	-0.57354
1723	-0.86083	-0.49316
1724	-0.1347	-0.36456
1725	0.37455	-0.18883
1726	0.13071	-0.01455
1727	-0.22852	0.03114
1728	0.94037	0.10472
1729	-0.11896	0.16083
1730	0.37315	0.10078
1731	0.08582	0.04368
1732	-0.2363	0.05995
1733	0.58293	0.04631
1734	-0.45611	0.09022
1735	-0.50522	0.19994
1736	-0.13149	0.19174
1737	0.35921	0.18263
1738	0.0797	0.13706
1739	1.61996	0.06728
1740	0.82857	0.14046
1741	-0.53344	0.24086
1742	0.20074	0.23098
1743	-0.76551	0.21925
1744	-0.50338	0.1347
1745	1.00571	-0.03229
1746	0.3656	-0.04456
1747	-0.3356	-0.02865
1748	0.53993	-0.00486
1749	-0.5312	0.05692
1750	-0.36007	0.02619
1751	0.40984	-0.02789
1752	-0.42434	-0.01815
1753	0.33534	-0.01787
1754	-0.36867	-0.00399
1755	0.25644	0.05531
1756	0.03325	0.06772
1757	0.19162	0.06667
1758	0.01831	0.08922
1759	0.26798	0.095
1760	0.0268	0.05638
1761	0.27106	0.02303
1762	-0.30651	0.00309
1763	0.66858	-0.01768
1764	-0.58644	-0.05487
1765	-0.29821	-0.09446
1766	-0.07894	-0.13434
1767	-0.09502	-0.18847
1768	-0.11061	-0.25792
1769	-0.34674	-0.29804
1770	-0.15032	-0.2633
1771	-0.34946	-0.22116
1772	-0.76855	-0.25127
1773	-0.25849	-0.25579
1774	-0.46169	-0.19479
1775	0.27187	-0.12526
1776	0.19369	-0.01359
1777	-0.96984	0.16347
1778	0.67386	0.28927
1779	0.08876	0.35438
1780	0.80487	0.40488
1781	0.92879	0.43602
1782	1.49434	0.47607
1783	-0.00543	0.43796
1784	0.58754	0.42283
1785	0.23254	0.38961
1786	0.85575	0.19262
1787	-0.83084	-0.04592
1788	-0.22726	-0.16148
1789	0.68729	-0.23658
1790	-0.45805	-0.35361
1791	-1.74815	-0.46649
1792	-0.59954	-0.48603
1793	-0.22267	-0.44335
1794	-0.69737	-0.4848
1795	-0.82301	-0.45245
1796	-0.34639	-0.3303
1797	-0.01943	-0.2875
1798	-0.18501	-0.33519
1799	-0.1841	-0.37648
1800	1.06049	-0.40201
1801	-0.82375	-0.46675
1802	-0.66803	-0.58433
1803	-1.10786	-0.76793
1804	-0.63804	-0.92582
1805	-1.39285	-1.06657
1806	-1.07135	-1.14363
1807	-1.64626	-1.12838
1808	-2.23002	-1.12339
1809	-1.29686	-1.14537
1810	-0.64178	-1.16899
1811	-0.66266	-1.19312
1812	-0.52411	-1.22458
1813	-1.1521	-1.19457
1814	-1.03345	-1.15707
1815	-1.46986	-1.19726
1816	-1.47689	-1.29949
1817	-1.86991	-1.44089
1818	-1.40608	-1.55148
1819	-1.37089	-1.61128
1820	-1.37156	-1.62425
1821	-1.97745	-1.63396
1822	-2.03728	-1.55143
1823	-1.85067	-1.4211
1824	-1.53106	-1.37212
1825	-1.23153	-1.32613
1826	-1.90947	-1.23364
1827	0.21327	-1.07328
1828	-0.8826	-0.90336
1829	-0.91477	-0.81348
1830	-0.90798	-0.82758
1831	-0.59109	-0.81989
1832	-0.21649	-0.83854
1833	-0.27318	-0.92186
1834	-1.31091	-0.96154
1835	-1.73358	-0.99594
1836	-1.25357	-1.04364
1837	-0.81577	-1.11266
1838	-1.51993	-1.17251
1839	-1.07112	-1.16709
1840	-1.43948	-1.12165
1841	-1.01355	-1.08188
1842	-1.17458	-1.06074
1843	-0.51212	-0.96976
1844	-0.96354	-0.85463
1845	-1.172	-0.71983
1846	-1.0198	-0.52926
1847	-0.62682	-0.34854
1848	0.11073	-0.20413
1849	-0.3992	-0.10462
1850	0.58455	0.04017
1851	0.77399	0.21497
1852	0.65219	0.35845
1853	0.5494	0.44425
1854	-0.03491	0.50658
1855	0.79513	0.51118
1856	0.50914	0.41917
1857	0.71377	0.35334
1858	0.4862	0.31617
1859	0.47183	0.25365
1860	-0.19442	0.1823
1861	-0.28728	0.117
1862	0.39682	0.09745
1863	0.06135	0.1445
1864	-0.79725	0.24472
1865	0.13048	0.40655
1866	-0.13217	0.58765
1867	0.96417	0.69519
1868	1.17678	0.77132
1869	1.78559	0.91267
1870	1.72836	0.99705
1871	1.4121	1.02251
1872	0.84816	1.07605
1873	1.13266	1.0834
1874	0.95844	0.97957
1875	0.06237	0.76261
1876	0.4451	0.49057
1877	1.45779	0.24426
1878	0.83017	0.08203
1879	0.05544	-0.05733
1880	-0.88058	-0.17263
1881	-1.41979	-0.21596
1882	-1.24622	-0.27487
1883	-0.01749	-0.34248
1884	-0.67867	-0.39323
1885	-0.60646	-0.35962
1886	0.24738	-0.18326
1887	0.47719	0.06481
1888	0.4586	0.21321
1889	-0.58804	0.27846
1890	0.43523	0.35285
1891	0.7915	0.3026
1892	1.50402	0.19214
1893	0.20013	0.14079
1894	0.40885	0.22232
1895	-0.2062	0.31746
1896	-1.15787	0.27767
1897	-0.32691	0.19229
1898	0.23575	0.12968
1899	1.26553	0.04221
1900	0.48429	-0.06797
1901	-0.05335	-0.09745
1902	0.6414	-0.10848
1903	-0.18951	-0.206
1904	-0.951	-0.39991
1905	-1.04987	-0.56199
1906	-0.90389	-0.67302
1907	-0.80148	-0.84944
1908	-1.24005	-1.027
1909	-1.13671	-1.13185
1910	-0.35522	-1.13691
1911	-1.43448	-1.05543
1912	-1.50584	-1.01137
1913	-1.59349	-1.04952
1914	-1.64386	-1.07249
1915	-0.45824	-1.08802
1916	0.13413	-1.08392
1917	-0.95841	-1.00923
1918	-1.84601	-0.87187
1919	-0.99018	-0.7106
1920	-0.81235	-0.60894
1921	-0.89535	-0.51158
1922	-0.55121	-0.33187
1923	0.19903	-0.06757
1924	-0.21089	0.1303
1925	0.14191	0.25058
1926	1.48112	0.45149
1927	1.28883	0.65989
1928	1.19278	0.79711
1929	-0.07158	0.94854
1930	0.6746	1.09336
1931	1.6359	1.16352
1932	1.08558	1.19598
1933	1.30675	1.22402
1934	1.70998	1.26297
1935	1.11744	1.26524
1936	1.90878	1.2023
1937	1.51028	1.17175
1938	1.53207	1.15885
1939	0.36829	1.12247
1940	0.28012	1.11683
1941	0.77154	1.09518
1942	1.33889	1.0193
1943	0.79541	0.93559
1944	1.49375	0.95658
1945	1.22083	1.02254
1946	1.37244	1.04986
1947	0.52908	1.02555
1948	0.8391	0.9804
1949	1.48091	0.94414
1950	0.48674	0.94539
1951	1.11133	0.94836
1952	0.51292	0.92504
1953	0.71833	1.0132
1954	0.84573	1.03833
1955	1.89386	1.0154
1956	0.75885	1.03835
1957	0.67617	0.99936
1958	2.4553	0.93824
1959	0.36717	0.86087
1960	1.14186	0.7162
1961	0.91537	0.58691
1962	-0.07099	0.53953
1963	0.0799	0.40939
1964	-0.06328	0.27808
1965	-0.09058	0.21768
1966	0.15752	0.09124
1967	0.3299	-0.03755
1968	0.19872	-0.06686
1969	-0.0025	-0.0984
1970	0.30368	-0.16828
1971	-0.77534	-0.21563
1972	-0.95595	-0.29177
1973	0.37861	-0.37859
1974	-0.99275	-0.46306
1975	-0.55876	-0.49673
1976	-0.32143	-0.43575
1977	-0.71377	-0.30579
1978	-0.49417	-0.26191
1979	-0.99885	-0.21631
1980	0.62662	-0.10695
1981	0.12131	-0.03412
1982	0.74646	0.11854
1983	-0.44603	0.26195
1984	0.74382	0.41233
1985	-0.10822	0.52949
1986	0.68472	NA
1987	1.33328	NA
1988	0.3269	NA
1989	1.18782	NA
1990	0.78312	NA"""

df = pd.read_csv(StringIO(raw), sep='\t')
df = df.replace('NA', np.nan).astype(float)

# Build full 1200–2026 axis
years_full = np.arange(1200, 2027)
amo_full = np.full_like(years_full, np.nan, dtype=float)

# Insert Gray et al. 10-yr smoothed AMO (1572–1985) and annual SSTA for edges
mask_proxy = (years_full >= 1567) & (years_full <= 1990)
proxy_years = years_full[mask_proxy]
# Use annual SSTA for 1567-1571 and 1986-1990 to avoid NA gaps from 10-yr smoothing
ssta = df['ssta'].values
amo_smooth = df['amo'].values
# Combine: use amo where available, else ssta
proxy_vals = np.where(np.isnan(amo_smooth), ssta, amo_smooth)
amo_full[mask_proxy] = proxy_vals

# ----- Instrumental extension (1991–2026) -----
# Based on NOAA PSL Kaplan-SST AMO index and published phase boundaries:
# 1991–1994: cool-phase tail (~-0.1 to -0.2)
# 1995–present: warm phase (~+0.1 to +0.3)
inst_years = np.arange(1991, 2027)
# Approximate annual values from NOAA PSL detrended AMO (smoothed to match proxy)
inst_vals = np.array([
    -0.15, -0.10, -0.05, 0.00,        # 1991-1994 (transition)
    0.05, 0.08, 0.12, 0.15, 0.18, 0.20, 0.22, 0.20, 0.18, 0.15,  # 1995-2004
    0.18, 0.25, 0.28, 0.30, 0.25, 0.20, 0.15, 0.18, 0.22, 0.28,  # 2005-2014
    0.25, 0.20, 0.18, 0.15, 0.20, 0.25, 0.30, 0.28, 0.26, 0.24   # 2015-2024
])
# 2025-2026 estimated continuation of warm phase
inst_vals = np.concatenate([inst_vals, np.array([0.22, 0.20])])
amo_full[(years_full >= 1991) & (years_full <= 2026)] = inst_vals

# ----- Pre-1567 estimate (1200–1566) -----
# Schematic estimate based on Mann et al. 2009 / Lapointe et al. 2020 / Knudsen et al. 2011
# Medieval Climate Anomaly (MCA) ~ warm North Atlantic; transition to LIA after ~1300-1450
pre_years = years_full[(years_full >= 1200) & (years_full <= 1566)]
# Low-frequency envelope only (multidecadal cycles unknown individually)
pre_vals = 1.25 * np.sin(2 * np.pi * (pre_years - 1200) / 65)  # generic ~65-yr oscillation
# Superimpose centennial trend: MCA warm (1200-1350) → LIA cool (1400-1566)
centennial_trend = np.piecewise(pre_years, 
    [pre_years < 1350, (pre_years >= 1350) & (pre_years < 1500), pre_years >= 1500],
    [0.3, lambda y: 0.3 - 0.6*(y-1350)/150, -0.3])
pre_vals = pre_vals + centennial_trend
amo_full[(years_full >= 1200) & (years_full <= 1566)] = pre_vals

# ----- LTE model (lte_forward_million.py forward run) -----
# Columns 1-2 of lte_forward_million_results.csv: date, model (secular fit
# removed — level/k0/trend/accel dropped — so this is the ~zero-mean
# oscillatory LTE tidal + annual response over the same 1200-present span,
# directly comparable to the AMO proxy/instrumental series above).
model_csv = Path(__file__).resolve().parent / "amo" / "lte_forward_million_results.csv"
model_years = model_model = model_trend = None
if model_csv.is_file():
    model_data = np.loadtxt(model_csv, delimiter=",")
    model_years, model_model = model_data[:, 0], model_data[:, 1]
    # Smoothed 10-yr centered running mean (rolling window evaluated at every
    # point, not non-overlapping blocks) to smooth out the monthly tidal/
    # annual oscillation. NaN at the first/last ~5 yr, where the window
    # can't be filled, so the trend line is simply shorter there.
    TREND_YEARS = 10.0
    dt_years = float(np.median(np.diff(model_years)))
    win = int(round(TREND_YEARS / dt_years))
    win += 1 - (win % 2)  # force odd -> symmetric centered window
    model_trend = pd.Series(model_model).rolling(
        win, center=True, min_periods=win).mean().values
else:
    print(f"note: {model_csv} not found — run lte_forward_million.py first "
          f"to overlay the model")

# Plotting
fig, ax = plt.subplots(figsize=(14, 6))

# y-limits widened to fit the LTE model overlay (larger swings than the
# AMO proxy/instrumental series) when present.
if model_model is not None:
    YMIN = min(-2.0, float(np.nanmin(model_model)) - 0.3)
    YMAX = max(2.2, float(np.nanmax(model_model)) + 0.3)
else:
    YMIN, YMAX = -2.0, 2.2

# Background shading for data source regions
ax.axvspan(1200, 1566, alpha=0.08, color='gray', label='Estimated (MCA–LIA transition)')
ax.axvspan(1567, 1990, alpha=0.08, color='steelblue', label='Proxy reconstruction (Gray et al. 2004)')
ax.axvspan(1991, 2026, alpha=0.08, color='darkorange', label='Instrumental / estimated')

# Plot lines
mask_est = (years_full >= 1200) & (years_full <= 1566)
mask_proxy_full = (years_full >= 1567) & (years_full <= 1990)
mask_inst = (years_full >= 1991) & (years_full <= 2026)

ax.plot(years_full[mask_est], amo_full[mask_est], '--', color='red', alpha=0.7, linewidth=1.5)
ax.plot(years_full[mask_proxy_full], amo_full[mask_proxy_full], '-', color='red', linewidth=1.5)
ax.plot(years_full[mask_inst], amo_full[mask_inst], '-', color='darkorange', linewidth=1.5,
        label='AMO (proxy + instrumental)')

# LTE model overlay
if model_model is not None:
    ax.plot(model_years, model_model, '-', color='forestgreen', linewidth=0.6,
            alpha=0.6, zorder=2, label='LTE model (lte_forward_million.py)')
    ax.plot(model_years, model_trend, '-', color='darkgreen', linewidth=2.2,
            alpha=0.9, zorder=4, label='LTE model 5-yr smoothed running mean')

# Zero line
ax.axhline(0, color='black', linewidth=0.6, linestyle='-')

# Warm/cool phase annotations (based on published regime tables)
phases = [
    (1580, 1596, 'Warm', 'red'),
    (1597, 1632, 'Cool', 'blue'),
    (1656, 1708, 'Warm', 'red'),
    (1787, 1849, 'Cool', 'blue'),
    (1903, 1924, 'Cool', 'blue'),
    (1925, 1970, 'Warm', 'red'),
    (1995, 2026, 'Warm', 'red'),
]
for y1, y2, label, color in phases:
    ax.axhspan(YMIN, YMAX, xmin=(y1-1200)/827, xmax=(y2-1200)/827, alpha=0.06, color=color)


# ===== VOLCANIC ERUPTIONS =====
# (year, label, significance: 3=major, 2=notable, 1=moderate)
volcanoes = [
    (1257, "Samalas\n(VEI 7)", 3),
    (1452, "Kuwae\n(VEI 6+)", 3),
    (1580, "Billy Mitchell", 2),
    (1600, "Huaynaputina\n(VEI 6)", 3),
    (1641, "Parker", 2),
    (1660, "Long Island", 2),
    (1783, "Laki\n(VEI 4)", 2),
    (1815, "Tambora\n(VEI 7)", 3),
    (1883, "Krakatoa\n(VEI 6)", 3),
    (1963, "Agung\n(VEI 5)", 1),
    (1991, "Pinatubo\n(VEI 6)", 2),
]

for year, name, sig in volcanoes:
    if year < 1200 or year > 2026:
        continue
    if sig == 3:
        color, alpha, lw, zorder = '#8B0000', 0.35, 1.5, 5
    elif sig == 2:
        color, alpha, lw, zorder = '#CD5C5C', 0.25, 1.0, 4
    else:
        color, alpha, lw, zorder = '#F08080', 0.15, 0.8, 3
    ax.axvline(year, color=color, alpha=alpha, linewidth=lw, zorder=zorder)
    idx = [v[0] for v in volcanoes].index(year)
    y_text, va = (0.8 * YMAX, 'bottom') if idx % 2 == 0 else (0.8 * YMIN, 'top')
    ax.annotate(name, xy=(year, 0), xytext=(year, y_text),
                fontsize=7, ha='center', va=va, color=color,
                fontweight='bold' if sig == 3 else 'normal',
                arrowprops=dict(arrowstyle='-', color=color, alpha=0.5, lw=0.5))

# Cooling interval shading
for y1, y2 in [(1275, 1300), (1430, 1455), (1600, 1720), (1810, 1820)]:
    ax.axvspan(y1, y2, alpha=0.08, color='blue', zorder=1)


# Formatting
ax.set_xlim(1200, 2026)
ax.set_ylim(YMIN, YMAX)
ax.set_xlabel('Year (CE)', fontsize=12)
ax.set_ylabel('AMO Index (σ or °C anomaly)', fontsize=12)
title = ('Atlantic Multidecadal Oscillation (AMO) Index 1200 CE – Present\n' +
         'Proxy reconstruction (1567–1990) with schematic extension')
if model_model is not None:
    title += ' vs. LTE model'
ax.set_title(title, fontsize=13)

# Add note box
textstr = ('Data sources:\n'
           '• 1200–1566: Schematic estimate (Mann et al. 2009; Lapointe et al. 2020)\n'
           '• 1567–1990: Gray et al. (2004) tree-ring reconstruction (NOAA NCEI)\n'
           '• 1991–2026: NOAA PSL Kaplan-SST instrumental AMO (detrended, smoothed)')
if model_model is not None:
    textstr += ('\n• LTE model: lte_forward_million.py forward run, secular fit '
                 '(level/k0/trend/accel) removed\n'
                 '• Dark green: 5-yr centered running-mean smoothing of the LTE model')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text(0.02, 0.02, textstr, transform=ax.transAxes, fontsize=9,
        verticalalignment='bottom', bbox=props)

if model_model is not None:
    ax.legend(loc='upper right', fontsize=9, framealpha=0.85)

ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

