# Mini Navigation & Positioning（迷你导航定位）

**从零开始、零依赖的 C 语言实现**，涵盖自主系统所需的导航、定位与授时算法。每个模块对应 MIT、Stanford 及其他顶尖大学的一门或多门课程，将教科书中的公式转化为可运行的 C 代码，实现理论与实践的桥接。

## 子模块总览

| 子模块 | 主题 | 参考课程 |
|--------|--------|-------------|
| [mini-geomagnetic-navigation](mini-geomagnetic-navigation/) | IGRF/WMM 球谐地磁场模型、磁图匹配（MAGCOM）、梯度导航、磁力计硬铁/软铁标定、卡尔曼滤波磁辅助惯导 | MIT 16.687, Stanford AA272C |
| [mini-gnss-gps](mini-gnss-gps/) | ECEF/大地坐标变换、伪距修正、Bancroft 直接解算、迭代最小二乘 PVT、Hatch 载波相位平滑、Klobuchar 电离层模型、Saastamoinen 对流层模型、DOP、C/A 码生成、多普勒、捕获搜索 | Stanford AA272C, MIT 16.687, Misra & Enge |
| [mini-indoor-positioning](mini-indoor-positioning/) | WiFi/BLE/磁指纹定位、行人航位推算（PDR）、IMU 惯性导航、ToF/TDoA/AOA 时间法定位、线性和扩展卡尔曼滤波传感器融合、误差指标（DOP、CEP、RMSE） | CMU 16-833, Stanford CS225A, MIT 6.882 |
| [mini-inertial-navigation](mini-inertial-navigation/) | 四元数/方向余弦阵/欧拉角姿态表示、IMU 标定（零偏、比例因子、非正交误差）、捷联力学编排方程、舒勒振荡、Allan 方差误差表征、GNSS/INS 松组合 | MIT 16.485, Stanford AA272C, Titterton & Weston |
| [mini-integrated-navigation](mini-integrated-navigation/) | GNSS 加权最小二乘定位解算、卫星星历位置计算、IMU 力学编排、松/紧/深组合 GNSS+INS 架构、卡尔曼滤波框架、旋转表示方法 | Stanford AA272C, MIT 16.687, Groves (2013) |
| [mini-slam-basics](mini-slam-basics/) | EKF-SLAM 地标初始化、FastSLAM（Rao-Blackwellized 粒子滤波）、图优化位姿图、数据关联（马氏距离、JCBB）、距离-方位传感器模型、SE(2)/SE(3) 李群位姿 | MIT 16.485, Stanford CS231A, Thrun et al. |
| [mini-timing-sync](mini-timing-sync/) | Allan 方差与频率稳定度、时钟模型（偏移/偏斜/漂移/老化）、NTP 客户端（RFC 5905）、PTP 引擎（IEEE 1588-2019）、鉴相/鉴频器、数字锁相环、时间传递（双向、共视、全视） | NIST, IEEE 1588, MIT 6.241J |
| [mini-uwb-localization](mini-uwb-localization/) | IEEE 802.15.4a 信道模型（LOS/NLOS 住宅/办公/工业）、Saleh-Valenzuela 多径模型、双向测距（TWR）、TOA 定位、NLOS 检测与抑制、粒子滤波跟踪、自包含线性代数 | IEEE 802.15.4a, Stanford CS444M, MIT 6.882 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有课程对齐说明与参考文献注释
- **实用演示程序** — GNSS 接收机、INS 力学编排、SLAM 机器人、NTP 时间服务器、UWB 定位器

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-geomagnetic-navigation
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-navigation-positioning/
├── mini-geomagnetic-navigation/   # IGRF/WMM 球谐地磁模型、MAGCOM、磁辅助卡尔曼滤波
├── mini-gnss-gps/                 # 伪距、Bancroft PVT、电离层/对流层、C/A 码、DOP
├── mini-indoor-positioning/       # WiFi/BLE 指纹定位、PDR、ToF/TDoA 传感器融合
├── mini-inertial-navigation/      # 四元数/方向余弦阵、捷联力学编排、IMU 标定
├── mini-integrated-navigation/    # GNSS+INS 耦合架构、卡尔曼滤波框架
├── mini-slam-basics/              # EKF-SLAM、FastSLAM、位姿图优化、数据关联
├── mini-timing-sync/              # NTP、PTP、Allan 方差、时钟模型、数字锁相环
└── mini-uwb-localization/         # IEEE 802.15.4a 信道、TWR 测距、NLOS 抑制
```

## 许可证

MIT
