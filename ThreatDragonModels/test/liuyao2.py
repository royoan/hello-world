import json

# ===================== 全局常量定义（统一所有重复配置）=====================
# 基础权重配置
BASE_WEIGHTS = {
    "月建": 1.0, "月建空亡": 0.5, "日辰": 0.8, "动爻": 0.6,
    "变爻": 0.4, "元神静": 0.3, "忌神静": -0.3, "仇神静": -0.2
}

# 四时旺相系数
SI_SHI = {
    "木": {"春": 0.2, "夏": -0.1, "秋": -0.2, "冬": 0.1},
    "火": {"春": 0.1, "夏": 0.2, "秋": -0.1, "冬": -0.2},
    "土": {"春": -0.1, "夏": 0.1, "秋": 0.1, "冬": -0.1},
    "金": {"春": -0.2, "夏": -0.1, "秋": 0.2, "冬": 0.1},
    "水": {"春": 0.1, "夏": -0.2, "秋": -0.1, "冬": 0.2}
}

# 地支五行对应
DI_ZHI_WU_XING = {
    "子": "水", "丑": "土", "寅": "木", "卯": "木",
    "辰": "土", "巳": "火", "午": "火", "未": "土",
    "申": "金", "酉": "金", "戌": "土", "亥": "水"
}

# 地支相冲对应
DI_ZHI_CHONG = {
    "子": "午", "午": "子", "丑": "未", "未": "丑",
    "寅": "申", "申": "寅", "卯": "酉", "酉": "卯",
    "辰": "戌", "戌": "辰", "巳": "亥", "亥": "巳"
}

# 地支相合对应
DI_ZHI_HE = {
    "子": "丑", "丑": "子", "寅": "亥", "亥": "寅",
    "卯": "戌", "戌": "卯", "辰": "酉", "酉": "辰",
    "巳": "申", "申": "巳", "午": "未", "未": "午"
}

# 旬信息字典（标准六爻规则）
XUN_INFO = {
    "甲子": {"include": ["子","丑","寅","卯","辰","巳","午","未","申","酉"], "empty": ["戌","亥"]},
    "甲戌": {"include": ["戌","亥","子","丑","寅","卯","辰","巳","午","未"], "empty": ["申","酉"]},
    "甲申": {"include": ["申","酉","戌","亥","子","丑","寅","卯","辰","巳"], "empty": ["午","未"]},
    "甲午": {"include": ["午","未","申","酉","戌","亥","子","丑","寅","卯"], "empty": ["辰","巳"]},
    "甲辰": {"include": ["辰","巳","午","未","申","酉","戌","亥","子","丑"], "empty": ["寅","卯"]},
    "甲寅": {"include": ["寅","卯","辰","巳","午","未","申","酉","戌","亥"], "empty": ["子","丑"]}
}

# 五行相生（全局唯一定义）
SHENG = {"木": "火", "火": "土", "土": "金", "金": "水", "水": "木"}

# 五行相克（全局唯一定义）
KE = {"木": "土", "土": "水", "水": "火", "火": "金", "金": "木"}

# ===================== 六爻股票预测类 =====================
class LiuYaoStockPredictor:
    def __init__(self, params):
        """
        初始化六爻股票预测器
        params格式：
        {
            "求测事项": "招商银行股票今日涨跌",
            "起卦时间": "丙午年 庚寅月 己酉日",
            "本卦": "天火同人",
            "动爻": [("九五", "申", "妻财", "酉")],  # (爻位, 动爻地支, 六亲, 变爻地支)
            "变卦": None,
            "用神": {"爻位": "九五", "干支": "壬申金", "六亲": "妻财"},
            "世爻": {"爻位": "九三", "干支": "己亥水", "六亲": "官鬼"},
            "应爻": {"爻位": "上九", "干支": "壬戌土", "六亲": "子孙"},
            "卦中爻": [
                {"爻位": "初九", "干支": "己卯木", "六亲": "父母", "空亡": True},
                {"爻位": "六二", "干支": "己丑土", "六亲": "子孙", "空亡": False},
                {"爻位": "九三", "干支": "己亥水", "六亲": "官鬼", "空亡": False},
                {"爻位": "九四", "干支": "壬午火", "六亲": "兄弟", "空亡": False},
                {"爻位": "九五", "干支": "壬申金", "六亲": "妻财", "空亡": False},
                {"爻位": "上九", "干支": "壬戌土", "六亲": "子孙", "空亡": False}
            ],
            "特殊说明": ""
        }
        """
        self.params = params
        self._parse_time()
        self._parse_xun_kong()
        self._parse_yong_shen()
        self._parse_other_yao()

    def _parse_time(self):
        """解析起卦时间"""
        time_parts = self.params["起卦时间"].split()
        self.year = time_parts[0][:-1]
        self.month = time_parts[1][:-1]
        self.day = time_parts[2][:-1]
        self.month_gan = self.month[0]
        self.month_zhi = self.month[1]
        self.day_gan = self.day[0]
        self.day_zhi = self.day[1]
        
        # 确定季节
        if self.month_zhi in ["寅", "卯", "辰"]:
            self.season = "春"
        elif self.month_zhi in ["巳", "午", "未"]:
            self.season = "夏"
        elif self.month_zhi in ["申", "酉", "戌"]:
            self.season = "秋"
        else:
            self.season = "冬"

    def _parse_xun_kong(self):
        """解析旬空：根据日支判断所属旬首及空亡（符合《增删卜易》规则）"""
        for xun_shou, info in XUN_INFO.items():
            if self.day_zhi in info["include"]:
                self.xun_shou = xun_shou
                self.xun_kong_zhi = info["empty"]
                break
        else:
            self.xun_shou = "甲子"
            self.xun_kong_zhi = ["戌","亥"]
        
        # 判断月建是否空亡
        self.month_kong = self.month_zhi in self.xun_kong_zhi

    def _parse_yong_shen(self):
        """解析用神"""
        self.yong_shen = self.params["用神"]
        self.ys_zhi = self.yong_shen["干支"][1]
        self.ys_wu_xing = DI_ZHI_WU_XING[self.ys_zhi]
        self.ys_kong = any(yao["爻位"] == self.yong_shen["爻位"] and yao["空亡"] 
                          for yao in self.params["卦中爻"])
        self.ys_yue_po = DI_ZHI_CHONG[self.month_zhi] == self.ys_zhi

    def _parse_other_yao(self):
        """解析其他爻（元神/忌神/仇神）"""
        self.yuan_shen = []  
        self.ji_shen = []   
        self.chou_shen = [] 
        
        for yao in self.params["卦中爻"]:
            if yao["爻位"] == self.yong_shen["爻位"]:
                continue
            yao_zhi = yao["干支"][1]
            yao_wu_xing = DI_ZHI_WU_XING[yao_zhi]
            
            if SHENG[yao_wu_xing] == self.ys_wu_xing:
                self.yuan_shen.append({"data": yao, "kong": yao["空亡"]})
            elif KE[yao_wu_xing] == self.ys_wu_xing:
                self.ji_shen.append({"data": yao, "kong": yao["空亡"]})
            else:
                for ys in self.yuan_shen:
                    ys_yao_zhi = ys["data"]["干支"][1]
                    ys_yao_wu_xing = DI_ZHI_WU_XING[ys_yao_zhi]
                    if KE[yao_wu_xing] == ys_yao_wu_xing:
                        for js in self.ji_shen:
                            js_yao_zhi = js["data"]["干支"][1]
                            js_yao_wu_xing = DI_ZHI_WU_XING[js_yao_zhi]
                            if SHENG[yao_wu_xing] == js_yao_wu_xing:
                                self.chou_shen.append({"data": yao, "kong": yao["空亡"]})
                        break

    def calc_yue_jian_force(self):
        """计算月建作用力量（分层版：临月建/月生/月克/月合/月破）"""
        YUE_JIAN_WEIGHTS = {
            "临月建": 1.2,    "月建生用神": 1.0, "月建合用神": 0.9,
            "同类帮扶": 0.8,   "无作用": 0.0,     "月建克用神": -1.0,
            "月破": -1.5       
        }
        
        month_wu_xing = DI_ZHI_WU_XING[self.month_zhi]
        
        # 优先级：月破 > 临月建 > 月合 > 五行生克
        if DI_ZHI_CHONG[self.month_zhi] == self.ys_zhi:
            base_force = YUE_JIAN_WEIGHTS["月破"]
        elif self.month_zhi == self.ys_zhi:
            base_force = YUE_JIAN_WEIGHTS["临月建"]
        elif DI_ZHI_HE[self.month_zhi] == self.ys_zhi:
            base_force = YUE_JIAN_WEIGHTS["月建合用神"]
        else:
            if SHENG[month_wu_xing] == self.ys_wu_xing:
                base_force = YUE_JIAN_WEIGHTS["月建生用神"]
            elif KE[month_wu_xing] == self.ys_wu_xing:
                base_force = YUE_JIAN_WEIGHTS["月建克用神"]
            elif month_wu_xing == self.ys_wu_xing:
                base_force = YUE_JIAN_WEIGHTS["同类帮扶"]
            else:
                base_force = YUE_JIAN_WEIGHTS["无作用"]
        
        # 修正系数
        yue_kong_coeff = 0.5 if self.month_kong else 1.0
        si_shi_coeff = SI_SHI[self.ys_wu_xing][self.season]
        final_force = base_force * yue_kong_coeff * (1 + si_shi_coeff)
        
        return round(final_force, 4)

    def calc_ri_chen_force(self):
        """计算日辰作用力量（分层版：临日/日合/日生/日克/日冲）"""
        RI_CHEN_WEIGHTS = {
            "临日建": 1.2,       "日合（静爻合起）": 1.1, "日合（动爻合绊）": 0.6,
            "日生用神": 1.0,     "同类帮扶（拱扶）": 0.8, "无作用": 0.0,
            "日克用神": -1.0,    "日破（休囚静爻）": -1.5, "暗动（旺相静爻）": 0.9,
        }
        
        day_wu_xing = DI_ZHI_WU_XING[self.day_zhi]
        
        # 优先级：日冲 > 临日建 > 日合 > 五行生克
        if DI_ZHI_CHONG[self.day_zhi] == self.ys_zhi:
            # 区分暗动/日破
            is_wang = (DI_ZHI_WU_XING[self.month_zhi] == self.ys_wu_xing or 
                       SHENG[DI_ZHI_WU_XING[self.month_zhi]] == self.ys_wu_xing or 
                       self.month_zhi == self.ys_zhi)
            is_xiu = (KE[DI_ZHI_WU_XING[self.month_zhi]] == self.ys_wu_xing or 
                      DI_ZHI_CHONG[self.month_zhi] == self.ys_zhi)
            ys_is_dong = self.yong_shen["爻位"] in [d[0] for d in self.params.get("动爻", [])]
            
            if ys_is_dong:
                base_force = RI_CHEN_WEIGHTS["日克用神"]
            elif is_wang:
                base_force = RI_CHEN_WEIGHTS["暗动（旺相静爻）"]
            else:
                base_force = RI_CHEN_WEIGHTS["日破（休囚静爻）"]
        elif self.day_zhi == self.ys_zhi:
            base_force = RI_CHEN_WEIGHTS["临日建"]
        elif DI_ZHI_HE[self.day_zhi] == self.ys_zhi:
            ys_is_dong = self.yong_shen["爻位"] in [d[0] for d in self.params.get("动爻", [])]
            base_force = RI_CHEN_WEIGHTS["日合（动爻合绊）"] if ys_is_dong else RI_CHEN_WEIGHTS["日合（静爻合起）"]
        else:
            if SHENG[day_wu_xing] == self.ys_wu_xing:
                base_force = RI_CHEN_WEIGHTS["日生用神"]
            elif KE[day_wu_xing] == self.ys_wu_xing:
                base_force = RI_CHEN_WEIGHTS["日克用神"]
            elif day_wu_xing == self.ys_wu_xing:
                base_force = RI_CHEN_WEIGHTS["同类帮扶（拱扶）"]
            else:
                base_force = RI_CHEN_WEIGHTS["无作用"]
        
        # 空亡修正（日冲空为冲实）
        kong_coeff = 0.5 if (self.ys_kong and not DI_ZHI_CHONG[self.day_zhi] == self.ys_zhi) else 1.0
        final_force = base_force * kong_coeff
        
        return round(final_force, 4)

    def calc_dong_bian_force(self):
        """计算动变爻作用力量（终极版：含三合局/三刑局/化进退/回头生克）"""
        # 动爻配置
        DONG_YAO_STATUS = {"临月建": 1.2, "得月生": 1.0, "月破": 0.0, "临日建": 1.2, "得日生": 1.0, "日破": 0.5, "动空": 0.5, "正常": 1.0}
        HUA_BIAN_TYPES = {"化进神": 1.3, "化退神": 0.7, "回头生": 1.5, "回头克": 0.3, "化月破": 0.0, "化日破": 0.5, "化空亡": 0.5, "无化变": 1.0}
        DONG_YAO_TO_YONG_SHEN = {"生用神": 1.0, "克用神": -1.0, "冲用神": -1.2, "合用神": 0.8, "无作用": 0.0}
        SAN_HE_JU = {"申子辰水局": (["申", "子", "辰"], "水", 2.5), "亥卯未木局": (["亥", "卯", "未"], "木", 2.5), "寅午戌火局": (["寅", "午", "戌"], "火", 2.5), "巳酉丑金局": (["巳", "酉", "丑"], "金", 2.5)}
        SAN_XING_JU = {"寅巳申持势之刑": (["寅", "巳", "申"], -2.0), "丑戌未无恩之刑": (["丑", "戌", "未"], -1.8), "子卯无礼之刑": (["子", "卯"], -1.5), "自刑": (["辰", "午", "酉", "亥"], -1.2)}
        
        # 初始化
        total_force = 0.0
        dong_yao_list = self.params.get("动爻", [])
        if not dong_yao_list:
            return 0.0
        
        # 辅助函数
        def is_xun_kong(zhi):
            return zhi in XUN_INFO[self.xun_shou]["empty"]
        
        def get_hua_bian_type(dong_zhi, bian_zhi):
            if not bian_zhi:
                return "无化变"
            # 化进退
            jin_shen = [("亥", "子"), ("寅", "卯"), ("巳", "午"), ("申", "酉"), ("丑", "辰"), ("辰", "未"), ("未", "戌"), ("戌", "丑")]
            tui_shen = [("子", "亥"), ("卯", "寅"), ("午", "巳"), ("酉", "申"), ("辰", "丑"), ("未", "辰"), ("戌", "未"), ("丑", "戌")]
            if (dong_zhi, bian_zhi) in jin_shen:
                return "化进神"
            elif (dong_zhi, bian_zhi) in tui_shen:
                return "化退神"
            # 回头生克
            dong_wx = DI_ZHI_WU_XING[dong_zhi]
            bian_wx = DI_ZHI_WU_XING[bian_zhi]
            if SHENG[bian_wx] == dong_wx:
                return "回头生"
            elif KE[bian_wx] == dong_wx:
                return "回头克"
            # 化空破
            if DI_ZHI_CHONG.get(self.month_zhi) == bian_zhi:
                return "化月破"
            elif DI_ZHI_CHONG.get(self.day_zhi) == bian_zhi:
                return "化日破"
            elif is_xun_kong(bian_zhi):
                return "化空亡"
            return "无化变"
        
        def get_dong_yao_status(dong_zhi):
            status = "正常"
            # 月建
            if dong_zhi == self.month_zhi:
                status = "临月建"
            elif DI_ZHI_CHONG.get(self.month_zhi) == dong_zhi:
                status = "月破"
            elif SHENG[DI_ZHI_WU_XING[self.month_zhi]] == DI_ZHI_WU_XING[dong_zhi]:
                status = "得月生"
            # 日辰覆盖
            if dong_zhi == self.day_zhi:
                status = "临日建"
            elif DI_ZHI_CHONG.get(self.day_zhi) == dong_zhi:
                status = "日破"
            elif SHENG[DI_ZHI_WU_XING[self.day_zhi]] == DI_ZHI_WU_XING[dong_zhi]:
                status = "得日生"
            # 旬空
            if is_xun_kong(dong_zhi):
                status = f"{status}_动空" if status != "正常" else "动空"
            return status
        
        def check_san_he_ju():
            """检查三合局"""
            all_zhi = set()
            dong_bian_zhi = set()
            for _, dong_zhi, _, bian_zhi in dong_yao_list:
                all_zhi.add(dong_zhi)
                dong_bian_zhi.add(dong_zhi)
                if bian_zhi:
                    all_zhi.add(bian_zhi)
                    dong_bian_zhi.add(bian_zhi)
            all_zhi.add(self.month_zhi)
            all_zhi.add(self.day_zhi)
            
            for ju_name, (zhi_list, ju_wx, ju_weight) in SAN_HE_JU.items():
                if set(zhi_list).issubset(all_zhi):
                    # 成局状态
                    has_empty = any(is_xun_kong(zhi) for zhi in zhi_list if zhi in dong_bian_zhi)
                    has_po = any(DI_ZHI_CHONG.get(self.month_zhi) == zhi or DI_ZHI_CHONG.get(self.day_zhi) == zhi for zhi in zhi_list if zhi in dong_bian_zhi)
                    ju_status = "实局" if not (has_empty or has_po) else "待局"
                    
                    # 作用方向
                    if ju_wx == self.ys_wu_xing or SHENG[ju_wx] == self.ys_wu_xing:
                        ju_direction = 1.0
                    elif KE[ju_wx] == self.ys_wu_xing:
                        ju_direction = -1.0
                    else:
                        ju_direction = 0.0
                    
                    final_weight = ju_weight * ju_direction
                    if ju_status == "待局":
                        final_weight *= 0.5
                    
                    self.params["san_he_info"] = {"局类型": ju_name, "成局状态": ju_status, "力量": final_weight}
                    return final_weight
            return 0.0
        
        def check_san_xing_ju():
            """检查三刑局"""
            dong_bian_zhi = set()
            for _, dong_zhi, _, bian_zhi in dong_yao_list:
                dong_bian_zhi.add(dong_zhi)
                if bian_zhi:
                    dong_bian_zhi.add(bian_zhi)
            
            for xing_name, (zhi_list, xing_weight) in SAN_XING_JU.items():
                if xing_name == "自刑":
                    for zhi in zhi_list:
                        if list(dong_bian_zhi).count(zhi) >= 2:
                            self.params["san_xing_info"] = {"刑类型": xing_name, "力量": xing_weight}
                            return xing_weight
                else:
                    xing_zhi_in = [zhi for zhi in zhi_list if zhi in dong_bian_zhi]
                    if len(xing_zhi_in) >= 2:
                        # 解刑判断
                        has_he = False
                        for zhi1 in xing_zhi_in:
                            for zhi2 in xing_zhi_in:
                                if zhi1 != zhi2 and DI_ZHI_HE.get(zhi1) == zhi2:
                                    has_he = True
                                    break
                            if has_he:
                                break
                        if not has_he:
                            final_weight = xing_weight * 1.2 if len(xing_zhi_in) == 3 else xing_weight
                            self.params["san_xing_info"] = {"刑类型": xing_name, "力量": final_weight}
                            return final_weight
            return 0.0
        
        def get_relation_type(dong_zhi, ys_zhi):
            """判断动爻与用神关系"""
            if DI_ZHI_CHONG.get(dong_zhi) == ys_zhi:
                return "冲用神"
            elif DI_ZHI_HE.get(dong_zhi) == ys_zhi:
                return "合用神"
            dong_wx = DI_ZHI_WU_XING[dong_zhi]
            ys_wx = DI_ZHI_WU_XING[ys_zhi]
            if SHENG[dong_wx] == ys_wx:
                return "生用神"
            elif KE[dong_wx] == ys_wx:
                return "克用神"
            return "无作用"
        
        # 核心计算
        # 1. 三合局（最高优先级）
        san_he_force = check_san_he_ju()
        total_force += san_he_force
        
        # 2. 三刑局（次高优先级）
        san_xing_force = check_san_xing_ju()
        total_force += san_xing_force
        
        # 3. 单动爻计算
        for _, dong_zhi, liu_qin, bian_zhi in dong_yao_list:
            # 动爻状态
            dong_status = get_dong_yao_status(dong_zhi)
            status_coeff = DONG_YAO_STATUS.get(dong_status.split("_")[0], 1.0)
            if "动空" in dong_status:
                status_coeff *= DONG_YAO_STATUS["动空"]
            
            # 化变类型
            hua_bian_type = get_hua_bian_type(dong_zhi, bian_zhi)
            hua_bian_coeff = HUA_BIAN_TYPES.get(hua_bian_type, 1.0)
            
            # 与用神关系
            relation_type = get_relation_type(dong_zhi, self.ys_zhi)
            relation_coeff = DONG_YAO_TO_YONG_SHEN.get(relation_type, 0.0)
            
            # 六亲权重
            liu_qin_weight = 1.2 if liu_qin in ["原神", "忌神"] else (0.8 if liu_qin == "仇神" else 0.5)
            
            # 单爻力量
            single_force = status_coeff * hua_bian_coeff * relation_coeff * liu_qin_weight
            total_force += single_force
        
        # 4. 多动爻衰减
        if len(dong_yao_list) > 1 and san_he_force == 0:
            multi_coeff = 1.0 - (len(dong_yao_list) - 1) * 0.1
            total_force *= max(multi_coeff, 0.5)
        
        return round(total_force, 4)

    def calc_jing_yao_force(self):
        """计算静爻生克力量"""
        yuan_shen_force = sum(BASE_WEIGHTS["元神静"] * (0.5 if ys["kong"] else 1.0) for ys in self.yuan_shen)
        ji_shen_force = sum(BASE_WEIGHTS["忌神静"] * (0.5 if js["kong"] else 1.0) for js in self.ji_shen)
        chou_shen_force = sum(BASE_WEIGHTS["仇神静"] * (0.5 if cs["kong"] else 1.0) for cs in self.chou_shen)
        return round(yuan_shen_force + ji_shen_force + chou_shen_force, 4)

    def calc_special_state_force(self):
        """计算特殊状态修正力量"""
        force = 0.0
        if self.ys_kong:
            force -= 0.5
        if self.ys_yue_po:
            force -= 0.3
        return round(force, 4)

    def calc_total_force(self):
        """计算总旺衰力量"""
        f_yue = self.calc_yue_jian_force()
        f_ri = self.calc_ri_chen_force()
        f_dong = self.calc_dong_bian_force()
        f_jing = self.calc_jing_yao_force()
        f_special = self.calc_special_state_force()
        
        total = f_yue + f_ri + f_dong + f_jing + f_special
        return {
            "总得分": round(total, 2),
            "月建得分": f_yue,
            "日辰得分": f_ri,
            "动变得分": f_dong,
            "静爻得分": f_jing,
            "特殊得分": f_special
        }

    def judge_ji_xiong(self):
        """判断吉凶"""
        total = self.calc_total_force()
        f_total = total["总得分"]
        
        # 世应系数
        shi_yao = self.params["世爻"]
        ying_yao = self.params["应爻"]
        shi_zhi = shi_yao["干支"][1]
        ying_zhi = ying_yao["干支"][1]
        shi_wu_xing = DI_ZHI_WU_XING[shi_zhi]
        ying_wu_xing = DI_ZHI_WU_XING[ying_zhi]
        
        if SHENG[shi_wu_xing] == ying_wu_xing or SHENG[ying_wu_xing] == shi_wu_xing:
            k_shi_ying = 1.2
        elif KE[shi_wu_xing] == ying_wu_xing or KE[ying_wu_xing] == shi_wu_xing:
            k_shi_ying = -1.2
        else:
            k_shi_ying = 1.0
        
        # 量能系数
        zi_sun_count = len([yao for yao in self.params["卦中爻"] if yao["六亲"] == "子孙"])
        zi_sun_kong = any(yao["六亲"] == "子孙" and yao["空亡"] for yao in self.params["卦中爻"])
        zi_sun_dong = any(yao["爻位"] in [d[0] for d in self.params.get("动爻", [])] and yao["六亲"] == "子孙" for yao in self.params["卦中爻"])
        
        if zi_sun_dong:
            k_liang_neng = 1.3
        elif zi_sun_count >= 2:
            k_liang_neng = 1.1
        elif zi_sun_kong:
            k_liang_neng = 0.8
        else:
            k_liang_neng = 1.0
        
        # 吉凶总分
        g = f_total * k_shi_ying * k_liang_neng
        
        # 判断结果
        if g >= 0.8:
            result = "大阳"
            desc = f"涨幅≥2%，总得分：{f_total}，吉凶系数：{round(g,2)}"
        elif g >= 0.3:
            result = "小阳"
            desc = f"涨幅0%-2%，总得分：{f_total}，吉凶系数：{round(g,2)}"
        elif g > -0.3:
            result = "十字星"
            desc = f"震荡±1%，总得分：{f_total}，吉凶系数：{round(g,2)}"
        elif g > -0.8:
            result = "小阴"
            desc = f"跌幅0%-2%，总得分：{f_total}，吉凶系数：{round(g,2)}"
        else:
            result = "大阴"
            desc = f"跌幅≥2%，总得分：{f_total}，吉凶系数：{round(g,2)}"
        
        return {
            "吉凶结果": result,
            "详细描述": desc,
            "世应系数": round(k_shi_ying, 2),
            "量能系数": round(k_liang_neng, 2),
            "吉凶总分": round(g, 2)
        }

    def get_stock_advice(self):
        """获取完整操作建议"""
        ji_xiong = self.judge_ji_xiong()
        total_force = self.calc_total_force()
        
        advice = {
            "核心建议": [],
            "风险提醒": []
        }
        
        # 核心建议
        if ji_xiong["吉凶结果"] in ["大阳", "小阳"]:
            advice["核心建议"].append(f"今日{'可重仓' if ji_xiong['吉凶总分'] >= 0.8 else '可轻仓'}介入，目标涨幅{ji_xiong['吉凶结果'] == '大阳' and '≥2%' or '0%-2%'}")
            advice["核心建议"].append("建议设置止盈位，逢高减仓")
        elif ji_xiong["吉凶结果"] == "十字星":
            advice["核心建议"].append("今日以观望为主，可小仓位高抛低吸，快进快出")
        else:
            advice["核心建议"].append("今日建议减仓或空仓观望，避免抄底")
        
        # 风险提醒
        if self.ys_yue_po:
            advice["风险提醒"].append("用神临月破，上涨动力有限，谨防冲高回落")
        if self.ys_kong:
            advice["风险提醒"].append("用神空亡，需等待出空填实后方可重仓操作")
        if "san_xing_info" in self.params:
            advice["风险提醒"].append(f"卦中出现{self.params['san_xing_info']['刑类型']}，易有突发风险，操作需谨慎")
        
        # 合局提示
        if "san_he_info" in self.params:
            he_info = self.params["san_he_info"]
            advice["核心建议"].append(f"卦中形成{he_info['局类型']}（{he_info['成局状态']}），力量{he_info['力量']}，需重点关注")
        
        return {
            "旺衰分析": total_force,
            "吉凶判断": ji_xiong,
            "操作建议": advice
        }

# ===================== 测试案例 =====================
if __name__ == "__main__":
    # 案例1：基础案例（无三合/三刑）
    params1 = {
        "求测事项": "招商银行股票今日涨跌",
        "起卦时间": "丙午年 庚寅月 己酉日",
        "本卦": "天火同人",
        "动爻": [("九五", "申", "妻财", "酉")],  # 申金动化酉金（进神）
        "变卦": None,
        "用神": {"爻位": "九五", "干支": "壬申金", "六亲": "妻财"},
        "世爻": {"爻位": "九三", "干支": "己亥水", "六亲": "官鬼"},
        "应爻": {"爻位": "上九", "干支": "壬戌土", "六亲": "子孙"},
        "卦中爻": [
            {"爻位": "初九", "干支": "己卯木", "六亲": "父母", "空亡": True},
            {"爻位": "六二", "干支": "己丑土", "六亲": "子孙", "空亡": False},
            {"爻位": "九三", "干支": "己亥水", "六亲": "官鬼", "空亡": False},
            {"爻位": "九四", "干支": "壬午火", "六亲": "兄弟", "空亡": False},
            {"爻位": "九五", "干支": "壬申金", "六亲": "妻财", "空亡": False},
            {"爻位": "上九", "干支": "壬戌土", "六亲": "子孙", "空亡": False}
        ],
        "特殊说明": "基础案例"
    }
    
    # 案例2：三合局案例（寅午戌火局生用神）
    params2 = {
        "求测事项": "贵州茅台股票今日涨跌",
        "起卦时间": "丙午年 庚寅月 庚午日",
        "本卦": "离为火",
        "动爻": [("六五", "戌", "官鬼", "未")],  # 戌土动
        "变卦": None,
        "用神": {"爻位": "九四", "干支": "戊午火", "六亲": "妻财"},
        "世爻": {"爻位": "九三", "干支": "戊辰土", "六亲": "父母"},
        "应爻": {"爻位": "上九", "干支": "己巳火", "六亲": "官鬼"},
        "卦中爻": [
            {"爻位": "初九", "干支": "己卯木", "六亲": "兄弟", "空亡": False},
            {"爻位": "六二", "干支": "己酉金", "六亲": "子孙", "空亡": False},
            {"爻位": "九三", "干支": "戊辰土", "六亲": "父母", "空亡": False},
            {"爻位": "九四", "干支": "戊午火", "六亲": "妻财", "空亡": False},
            {"爻位": "六五", "干支": "戊戌土", "六亲": "官鬼", "空亡": False},
            {"爻位": "上九", "干支": "己巳火", "六亲": "官鬼", "空亡": False}
        ],
        "特殊说明": "寅午戌火局生用神"
    }
    
    # 案例3：三刑局案例（寅巳申持势之刑克用神）
    params3 = {
        "求测事项": "宁德时代股票今日涨跌",
        "起卦时间": "丙午年 庚寅月 壬申日",
        "本卦": "雷水解",
        "动爻": [("九二", "寅", "父母", "卯"), ("九五", "巳", "官鬼", "午")],  # 寅+巳+申（日建）三刑
        "变卦": None,
        "用神": {"爻位": "六三", "干支": "丙申金", "六亲": "妻财"},
        "世爻": {"爻位": "初九", "干支": "庚子水", "六亲": "子孙"},
        "应爻": {"爻位": "上六", "干支": "丁未土", "六亲": "兄弟"},
        "卦中爻": [
            {"爻位": "初九", "干支": "庚子水", "六亲": "子孙", "空亡": False},
            {"爻位": "九二", "干支": "庚寅木", "六亲": "父母", "空亡": False},
            {"爻位": "六三", "干支": "丙申金", "六亲": "妻财", "空亡": False},
            {"爻位": "九四", "干支": "辛亥水", "六亲": "子孙", "空亡": False},
            {"爻位": "九五", "干支": "癸巳火", "六亲": "官鬼", "空亡": False},
            {"爻位": "上六", "干支": "丁未土", "六亲": "兄弟", "空亡": False}
        ],
        "特殊说明": "寅巳申三刑克用神"
    }

    params4 = {
        "求测事项": "半导体ETF-编码159813 明日涨跌",
        "起卦时间": "丙午年 庚寅月 庚戌日",
        "本卦": "乾为天",
        "动爻": [("九三", "辰", "父母", "丑")],
        "变卦": "天泽履",
        "用神": {"爻位": "九二", "干支": "甲寅木", "六亲": "妻财"},
        "世爻": {"爻位": "上九", "干支": "壬戌土", "六亲": "父母"},
        "应爻": {"爻位": "九三", "干支": "壬辰土", "六亲": "父母"},
        "卦中爻": [
            {"爻位": "初九", "干支": "己子水", "六亲": "子孙", "空亡": False},
            {"爻位": "九二", "干支": "己寅木", "六亲": "妻财", "空亡": True},
            {"爻位": "九三", "干支": "己辰土", "六亲": "父母", "空亡": False},
            {"爻位": "九四", "干支": "壬午火", "六亲": "官鬼", "空亡": False},
            {"爻位": "九五", "干支": "壬申金", "六亲": "兄弟", "空亡": False},
            {"爻位": "上九", "干支": "壬戌土", "六亲": "父母", "空亡": False}
        ],
        "特殊说明": "无"
    }
    
    # 运行测试
    print("="*60)
    print("测试案例1：基础案例")
    predictor1 = LiuYaoStockPredictor(params1)
    result1 = predictor1.get_stock_advice()
    print(f"旺衰分析：{result1['旺衰分析']}")
    print(f"吉凶判断：{result1['吉凶判断']['吉凶结果']}")
    print(f"操作建议：{result1['操作建议']['核心建议']}")
    
    print("\n" + "="*60)
    print("测试案例2：三合局案例")
    predictor2 = LiuYaoStockPredictor(params2)
    result2 = predictor2.get_stock_advice()
    print(f"旺衰分析：{result2['旺衰分析']}")
    print(f"吉凶判断：{result2['吉凶判断']['吉凶结果']}")
    print(f"操作建议：{result2['操作建议']['核心建议']}")
    
    print("\n" + "="*60)
    print("测试案例3：三刑局案例")
    predictor3 = LiuYaoStockPredictor(params3)
    result3 = predictor3.get_stock_advice()
    print(f"旺衰分析：{result3['旺衰分析']}")
    print(f"吉凶判断：{result3['吉凶判断']['吉凶结果']}")
    print(f"操作建议：{result3['操作建议']['核心建议']}")
    print(f"风险提醒：{result3['操作建议']['风险提醒']}")
    
    print("\n" + "="*60)
    print("测试案例4：半导体ETF-编码159813 明日涨跌")
    predictor4 = LiuYaoStockPredictor(params4)
    result4 = predictor4.get_stock_advice()
    print(f"旺衰分析：{result4['旺衰分析']}")
    print(f"吉凶判断：{result4['吉凶判断']['吉凶结果']}")
    print(f"操作建议：{result4['操作建议']['核心建议']}")
    print(f"风险提醒：{result4['操作建议']['风险提醒']}")
